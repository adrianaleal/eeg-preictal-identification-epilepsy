function [features, feature_names, comp_times, comp_times_names] = ...
    univariate_nonlinear_features(input_signal, fs, plotFigure)


% Input:
% - input_signal: signal to analyse
% - fs: sampling frequency
% - plotFigure: flag to plot figure

% Output:
% - features: univariate nonlinear features
% - feature_names: names of features
% - comp_times: computation time for each group of features
% - comp_times_names: name of each group of features



%% Fractal Dimension
tic
% Higuchi Fractal Dimension (Acharya2013)
% to choose Kmax (Kawe2019, Ruiz-Padial2018):

% Sharma2017: kmax = 128, 23.6 s EEG window, fs = 173.610 Hz (4097 samples)
% Kawe2019: kmax = 30, 1 s EEG window, fs = 256 Hz (256 samples)
% Ruiz-Padial2018 : kmax = 40, 21 s EEG window, fs = 500 Hz (256 samples)

Kmax = 100;
HFD = Higuchi_FD(input_signal, Kmax);
% Ruiz-Padial2018: "The range of values for HFD lie between 1 and 2
% 1 for a simple curve and 2 for a randomly distributed curve that nearly
% fills the Euclidean 2D space."

ct_Higuchi_FD = toc;


%% Detrended fluctuation analysis
% it is related to both autocorrelation and fractal dimension

% alpha_range = [10 32; 32 length(input_signal)/10, 10 length(input_signal)/10; ];

alpha_range = [10 32; 32 length(input_signal)/10];
% with this alpha_range we can see that the segments with small
% fluctuations have a random walk like structure whereas segments with
% large fluctuations have a noise like structure.


tic
% (2) Run a monofractal DFA in order to check if input_signal is a noise
% like time series
[~, ~, F, alpha_vec] = ...
    monofractal_detrended_fluctuation_analysis(input_signal, alpha_range, ...
    plotFigure);
DFA_alpha1 = alpha_vec(1);
DFA_alpha2 = alpha_vec(2);
ct_DFA = toc;

% (2) Run a multifractal DFA
tic

convert2RandomWalk = 1;
if alpha_vec(1)>=1.2 && alpha_vec(1)<=1.8
    disp('*******************************************************************')
    disp('*************** DETRENDED FLUCTUATION ANALYSIS ********************')
    disp('*******************************************************************')
    disp('The time series is random walk like')
    convert2RandomWalk = 0;
end
% if alpha_vec(1)>=0.2 && alpha_vec(1)<0.8
% disp('The time series is noise like')

[Ht,Htbin,Ph,Dh,MFDFA_ms_peak_abcissa, MFDFA_ms_width, MFDFA_ms_symmetry] = ...
    multifractal_detrended_fluctuation_analysis(...
    input_signal, alpha_range, convert2RandomWalk, plotFigure);
% feature that will be used from multifractal DFA: abcissa of the
% multifractal spectrum peak
ct_MFDFA = toc;


%% Multifractal 1-D wavelet leader (WL) estimates
tic
[Dh,h] = dwtleader(input_signal, 'db4');
% The wavelet leaders technique works best for data with 8000 or more samples.
% Serrano2009: "Wavelet Leaders exploits the wavelet self-similarity
% structures combined with the Multiresolution Analysis scheme."

Dh = flip(Dh);
h = flip(h);
hmin = h(1);
Dmin = Dh(1);
hmax = h(end);
Dmax = Dh(end);

ind = find(Dh==max(Dh));
hq0 = h(ind);
hq02hmin = hq0-hmin;
hmax2hq0 = hmax-hq0;

hmax2hmin = hmax-hmin;
Dmax2Dmin = Dmax-Dmin;

Dq02Dmin = Dh(ind)-Dmin;
Dq02Dmax = Dh(ind)-Dmax;
ct_WLMF = toc;

if plotFigure
    figure()
    subplot(211)
    time = 1/fs:1/fs:length(input_signal)/fs;
    plot(time, input_signal)
    axis tight
    ylabel('uV')
    xlabel({'Time (s)', '(a)'})
    subplot(212)
    hp = plot(h,Dh,'b-o');
    hp(1).MarkerFaceColor = 'b';
    grid on;
    xlabel('h'); ylabel('D(h)');
    title('Multifractal Spectrum');
    text(hmin-0.05,Dmin-0.04,'$[h_{min}, D_{min}]$', 'Interpreter', 'latex')
    text(hmax-0.05,Dmax-0.04,'$[h_{max}, D_{max}]$', 'Interpreter', 'latex')
    text(hq0-0.04,Dh(ind)+0.04,'$[h_{0}, D_{0}]$', 'Interpreter', 'latex')
    axis([-0.1 0.6 0 1.1])
%     name2save = 'pat_402_seiz_1_segment_1_count_win_segment_1';
%     save([name2save '.mat'], 'input_signal', 'fs', 'Dh', 'h')
    % export_fig(gcf, '-dpdf', fullfile(cd, 'FiguresPaper', 'dwtleader_pat_402_seiz1.pdf'), '-painters', '-transparent')
end

WL_MF = [hmin; Dmin; hmax; Dmax; hq0; hmax2hmin; Dmax2Dmin; ...
    Dq02Dmin; Dq02Dmax; hq02hmin; hmax2hq0];

wl_mf_feat_names = {'WL_MF_hmin'; 'WL_MF_Dmin'; 'WL_MF_hmax'; ...
    'WL_MF_Dmax'; 'WL_MF_hq0'; 'WL_MF_hmax2hmin'; 'WL_MF_Dmax2Dmin'; ...
    'WL_MF_Dq02Dmin'; 'WL_MF_Dq02Dmax'; 'WL_MF_hq02hmin'; 'WL_MF_hmax2hq0'};


%% normalize signal
tic
input_signal_norm = (input_signal-mean(input_signal))/std(input_signal); % Varsavsky2011


%% Phase Space Reconstruction
% (1) Time delay was estimated using the first local minimum of the average
% mutual information (using the default value of 'HistogramBins' = 10)
% The time delay was estimated using the default 'MaxLag' = 10
% (2) Embedding dimension was estimated using the False Nearest Neighbor
% algorithm (using the default value of 'PercentFalseNeighbors' = 0.1)
% The embedding dimension was estimated using the default 'MaxDim' = 5

[attractorPSR, tauPSR, eDimPSR] = phaseSpaceReconstruction(input_signal_norm);
ct_PSR = toc;

% if plotFigure
%     figure()
%     phaseSpaceReconstruction(input_signal_norm, tauPSR, eDimPSR);
% end


%% Sample and approximate entropy (Acharya2012, Acharya2015, Acharya2018)
tic
m = 2;
r_tolerance = 0.2;
SampEn = sample_entropy(input_signal_norm, m, r_tolerance);

% The default value of Radius is, 0.2*variance(X), if X has a
% single column.
% m default value is 2
ApEn = approximateEntropy(input_signal_norm);

ct_ApEnSampEn = toc;


%% Largest Lyapunov Exponent (BouAssi2017, Acharya2018)
tic
expansion_range = [1 20];
LLE = lyapunovExponent(input_signal_norm, fs, 'Dimension', eDimPSR, ...
    'Lag', tauPSR, 'ExpansionRange', expansion_range);
ct_LLE = ct_PSR + toc;


%% Correlation Dimension

% estimates the correlation dimension of the uniformly sampled
% time-domain signal

tic
min_radius = 0.1;
[CD, ~, ~] = correlationDimension(input_signal_norm, ...
    'Dimension', eDimPSR, 'Lag', tauPSR, 'MinRadius', min_radius);
ct_CD = ct_PSR + toc;

% Ma2018: "The values of the CD range between zero and the value of
% embedding dimension and can be used to quantify the complex dynamics of
% the brain activity [19]. In a chaotic system, the CD usually shows a
% non-integer value larger than 1, indicating an increased complexity of
% system dimensionality."


%%
if plotFigure
    
    lyapunovExponent(input_signal_norm, fs, 'Dimension', eDimPSR, ...
        'Lag', tauPSR, 'ExpansionRange', expansion_range);
    % MinSeparation is the threshold value used to find the nearest
    % neighbor i* for a point i to estimate the largest Lyapunov exponent.
    % 'MinSeparation', ceil(fs/max(meanfreq(X,fs)))
    f1 = gcf;
    Np = 50;
    correlationDimension(input_signal_norm,'Dimension', eDimPSR, 'Lag', tauPSR, ...
        'MinRadius', min_radius, 'NumPoints', Np);
    f2 = gcf;
    
    figure()
    set(gcf,'units','normalized','outerposition',[0 0 0.6 1])
    subplot(221)
    if eDimPSR==2
        plot(attractorPSR(:,1), attractorPSR(:,2))
        xlabel('x(k)')
        ylabel(['x(k-' num2str(tauPSR) ')'])
    elseif eDimPSR==3
        % A three-dimensional embedding uses a delay of 2*tau for the third
        % dimension.
        plot3(attractorPSR(:,1), attractorPSR(:,2), attractorPSR(:,3))
        xlabel('x(k)')
        ylabel('x(k-\tau)')
        zlabel('x(k-2\tau)')
    end
    box on
    
    s2 = subplot(222);
    copy_figure(f1, s2)
    
    s3 = subplot(223);
    copy_figure(f2, s3)
    
    % saveas(gcf, fullfile(cd, 'featureExtractionImages', ...
    %     [pat_seiz_name '_largLyapunovExp']))
    % close all
end


%% Recurrence Quantification Analysis
tic
[rqa_stat, ~, ~] = recurrenceAnalysis(input_signal_norm, tauPSR, ...
    eDimPSR, plotFigure, attractorPSR);
ct_RQA = ct_PSR + toc;

rqa_feat_names = {'RQA_REC'; 'RQA_DET'; 'RQA_Lmax'; 'RQA_L'; 'RQA_ENT'; ...
    'RQA_LAM'; 'RQA_TT'};


%% Simplicity measure
% Source: Nigam & Priemer, "Accessing heart dynamics to estimate durations
% of heart sounds", in Physiological Measurement 2005

% Also in: AN ADAPTIVE APPROACH TO ABNORMAL HEART SOUND SEGMENTATION
% Also in: Wavelet Transform And Simplicity Based Heart Murmur Segmentation

tic
emb_matrix = (attractorPSR/sqrt(length(attractorPSR)))'*(attractorPSR/sqrt(length(attractorPSR)));
eigvalues = sort(eig(emb_matrix),'descend');
eigvalues = eigvalues./sum(eigvalues);
H_simpl = -sum(eigvalues.*log10(eigvalues));
simplicity = 1/(2^H_simpl);
ct_simplicity = toc;


%% Store names and values of features

features = [HFD; DFA_alpha1; DFA_alpha2; MFDFA_ms_peak_abcissa; ...
    MFDFA_ms_width; MFDFA_ms_symmetry; WL_MF; ApEn; SampEn; tauPSR; ...
    eDimPSR; LLE; CD; rqa_stat'; simplicity];

feature_names = [{'Higuchi_FD'; 'DFA_alpha1'; 'DFA_alpha2'; ...
    'MFDFA_ms_peak'; 'MFDFA_ms_width'; 'MFDFA_ms_symmetry'}; ...
    wl_mf_feat_names; ...
    {'ApEn'; 'SampEn'; 'tauPSR'; 'eDimPSR'; 'LLE'; 'CD'}; ...
    rqa_feat_names; {'Simplicity'}];

comp_times = [ct_Higuchi_FD; ct_DFA; ct_MFDFA; ct_WLMF; ct_ApEnSampEn; ...
    ct_LLE; ct_CD; ct_RQA; ct_simplicity];

comp_times_names = {'ct_Higuchi_FD'; 'ct_DFA'; 'ct_MFDFA'; 'ct_ApEnSampEn'; ...
    'ct_LLE'; 'ct_CD'; 'ct_RQA'; 'ct_simplicity'};

end