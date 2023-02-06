function [features, feature_names, comp_time_vector, comp_time_names] = ...
    multivariate_features(win_all_chan_data, fs, eegChannelsName)


% Input:
% - win_all_chan_data: matrix of [n_chans x n_samples] containing the data
%   for all channels (n_chan) in a 5 second window (n_samples)
% - fs: sampling frequency
% - eegChannelsName: cell array containing the EEG channel names

% Output:
% - features: multivariate features
% - feature_names: names of features
% - comp_times: computation time for each group of features
% - comp_times_names: name of each group of features


% READ:
% (1) Baboukani2019: A novel multivariate phase synchrony measure: Application
% to multichannel newborn EEG analysis
% (2) Kida2016:


%%

tic


n_chans = numel(eegChannelsName);

electrode_pairs = nchoosek(1:n_chans,2);
n_electrode_pairs = size(electrode_pairs,1);
T = size(win_all_chan_data,2);
maxlag = T/4;

% initialize connectivity matrices:

n_measures = 2000; % randomly defined

% initialize connectivity matrices
undirect_measures_matrix = ones(n_chans, n_chans, n_measures);
% with exception of PLI, all measures take value 1 when analysing each
% channel with itself (elect1 vs elect1)

% initialize connectivity matrix
undirect_measures_names = cell(100,1);

% initialize connectivity vector
undirect_measures_vector = zeros(n_electrode_pairs, n_measures);

% initialize vector to save computation time to obtain the connectivity
% vector for each connectivity measure
comp_time_vector = zeros(n_electrode_pairs, n_measures);


% define the frequency bands
delta_band = [0.5, 4];
theta_band = [4, 8];
alpha_band = [8, 13];
beta_band = [13, 30];
gamma1_band = [30, 47];
gamma2_band = [53, 75];
gamma3_band = [75, 97];
gamma4_band = [103, 127];
freq_band_names = {'delta', 'theta', 'alpha', 'beta', 'gamma1', ...
    'gamma2', 'gamma3', 'gamma4'};
freq_bands = [delta_band; theta_band; alpha_band; beta_band; ...
    gamma1_band; gamma2_band; gamma3_band; gamma4_band];

n_bands = size(freq_bands, 1);


count_matrices = 0;

plotPhaseAnalysis = 0;

init_time = toc;

% compute the networks which are represented by their connectivity
% (adjacency) matrices:


for ee = 1:n_electrode_pairs
    
    tic
    count_measures = 0;
    
    elect1 = win_all_chan_data(electrode_pairs(ee,1),:);
    elect2 = win_all_chan_data(electrode_pairs(ee,2),:);
    get_elec_time = toc;
    
    % elect1 = (elect1-mean(elect1))/std(elect1);
    % elect2 = (elect2-mean(elect2))/std(elect2);
    
    for ff = 1:n_bands
        
        tic
        fc1 = freq_bands(ff,1);
        fc2 = freq_bands(ff,2);
        [z_ch1, amplitude_ch1, phi_ch1] = getAnalyticSignal(elect1, fs, ...
            fc1, fc2);
        [z_ch2, amplitude_ch2, phi_ch2] = getAnalyticSignal(elect2, fs, ...
            fc1, fc2);
        get_analytic_time = toc;
        
        %% get Spearman's correlation coefficient for instantaneous correlations
        tic
        RHO = corr(amplitude_ch1', amplitude_ch2', 'Type', 'Spearman');
        % nonparametric measure
        count_matrices = count_matrices+1;
        count_measures = count_measures+1;
        undirect_measures_matrix(electrode_pairs(ee,1), electrode_pairs(ee,2), count_measures) = RHO;
        undirect_measures_matrix(electrode_pairs(ee,2), electrode_pairs(ee,1), count_measures) = RHO;
        undirect_measures_vector(ee, count_measures) = RHO;
        undirect_measures_names(count_measures) = {['corr_coef_Spearman_' freq_band_names{ff}]};
        comp_time_vector(ee, count_measures) = toc+get_analytic_time+get_elec_time+init_time;
        
        %% get normalized cross-correlation (Mormann2003a, Mormann2005,
        % Mirowski2009)
        % cross-correlation analysis reveals whether peak connectivity is
        % observed when one time series is temporally shifted relative to
        % the other.(Cohen2014)
        tic
        max_r = max(xcorr(amplitude_ch1', amplitude_ch2', 'normalized'));
        
        %     [r, lags] = xcorr(elect1,elect1,maxlag);
        %     % unbised normalization:
        %     r = r./(T-abs(lags));
        %     % normalized normalization:
        %     cxx0 = sum(abs(elect1).^2);
        %     cyy0 = sum(abs(elect1).^2);
        %     scaleCoeffCross = sqrt(cxx0*cyy0);
        %     r = r./scaleCoeffCross;
        %     max_r2 = max(abs(r));
        
        count_matrices = count_matrices+1;
        count_measures = count_measures+1;
        undirect_measures_matrix(electrode_pairs(ee,1), electrode_pairs(ee,2), count_measures) = max_r;
        undirect_measures_matrix(electrode_pairs(ee,2), electrode_pairs(ee,1), count_measures) = max_r;
        undirect_measures_vector(ee, count_measures) = max_r;
        undirect_measures_names(count_measures) = {['corr_max_xcorr_' freq_band_names{ff}]};
        comp_time_vector(ee, count_measures) = toc+get_analytic_time+get_elec_time+init_time;
        
        %% get mutual information
%         plotFigure = 0;
%         ami = average_mutual_information([amplitude_ch1', amplitude_ch2'], plotFigure);
%         % [ MI ] = MI_KNN_cont_cont( amplitude_ch1', amplitude_ch2', 5 )
%         
%         count_matrices = count_matrices+1;
%         count_measures = count_measures+1;
%         undirect_measures_matrix(electrode_pairs(ee,1), electrode_pairs(ee,2), count_measures) = ami;
%         undirect_measures_matrix(electrode_pairs(ee,2), electrode_pairs(ee,1), count_measures) = ami;
%         undirect_measures_vector(ee, count_measures) = ami;
%         undirect_measures_names(count_measures) = {['Mutual_information_' freq_band_names{ff}]};
        
        %% get spectral coherence (magnitude-squared coherence) estimate
        %     % The magnitude-squared coherence enables the identification of
        %     % significant frequency-domain correlation between two time series.
        %     % The magnitude-squared coherence estimate is a function of frequency
        %     % with values between 0 and 1. These values indicate how well chan1
        %     % corresponds to chan2 at each frequency. The magnitude-squared
        %     % coherence is a function of the power spectral densities, Pxx(f) and
        %     % Pyy(f), and the cross power spectral density, Pxy(f), of x (elect1)
        %     % and y (elect2):
        %
        %     win = 1*fs;
        %     noverlap = [];
        %     nfft = [];
        %     [Cxy,F] = mscohere(elect1, elect2, win, noverlap, nfft, fs);
        %     % mscohere uses a Hamming window such that x and y are divided into
        %     % eight segments with noverlap overlapping samples.
        %     % with noverlap overlapping samples. window = length(elect1)/8 = 160
        %     % noverlap default: 50% = max(256,2^log2(win))
        %     % nfft default: max(256,2^log2(win))
        %     % nfft = max(256,2^nextpow2(160)) = 256
        %     % frequency resolution = 128 Hz / 256 FFT bins  = 0.5 Hz
        %
        %     % Cxy is similar to ISPC, but the phase values are weighted by power
        %     % values.
        %     % Cxy = 1 --> complete coherence
        %     % Cxy = 0 --> complete independence
        %     % despite being normalized by the total power in the denominator,
        %     % spectral coherence results can still be ingluenced by robust power
        %     % changes (Cohen2014, chapter 16)
        %
        %     figure()
        %     plot(F,Cxy)
        %     title('Coherence Estimate via Welch')
        %     ylabel('Magnitude-Squared Coherence')
        %     xlabel('Frequency (Hz)')
        %     grid
        %     axis tight
        %
        %     % autospectral power density for chan1
        %     [Pxx,~] = cpsd(elect1, elect1, [], [], [], fs);
        %     % autospectral power density for chan2
        %     [Pyy,~] = cpsd(elect2, elect2, [], [], [], fs);
        %     % cross power spectral density between electrodes
        %     [Pxy,F] = cpsd(elect1, elect2, [], [], [], fs);
        %
        %     Pxy(Cxy < 0.2) = 0;
        %
        %     figure()
        %     plot(F,angle(Pxy)/pi)
        %     title('Cross Spectrum Phase')
        %     xlabel('Frequency (Hz)')
        %     ylabel('Lag (\times\pi rad)')
        %     grid
        %     axis tight
        %
        %     Cxy2 = Pxy.*conj(Pxy)./(Pxx.*Pyy);
        %     delta_phi_cs = atan2(imag(Pxy), real(Pxy));
        %     delta_phi_cs2 = unwrap(angle(Pxy));
        %     delta_phi_cs3 = angle(Pxy)/pi;
        %
        %
        %     figure()
        %     plot(F, delta_phi_cs)
        %
        %
        %     ISPC_COH = abs(mean(exp(1i*angle(Pxy))));
        %
        %     % take imaginary part of signal only
        %     Pxyi = imag(Pxy);
        %
        %     % phase-lag index
        %     PLI_COH = abs(mean(sign(Pxyi)));
        
        
        %% get phase measures:
        tic
        % (3) shift phase angles:
        phase_ch1 = unwrap(phi_ch1);
        phase_ch2 = unwrap(phi_ch2);
        
        % (4) circular correlation between ch1 and ch2 (Baboukani2019):
        cir = circ_corrcc(phase_ch1, phase_ch2);
        % if cir = 1 --> ch1 = ch2+const (mod 2pi)
        % if cir = -1 --> ch1+ch2 = const (mod 2pi)
        % read Topics in Circular Statistics page 17
        
        count_matrices = count_matrices+1;
        count_measures = count_measures+1;
        undirect_measures_matrix(electrode_pairs(ee,1), electrode_pairs(ee,2), count_measures) = cir;
        undirect_measures_matrix(electrode_pairs(ee,2), electrode_pairs(ee,1), count_measures) = cir;
        undirect_measures_vector(ee, count_measures) = cir;
        undirect_measures_names(count_measures) = {['Circular_correlation_' freq_band_names{ff}]};
        comp_time_vector(ee, count_measures) = toc+get_analytic_time+get_elec_time+init_time;
 
        if plotPhaseAnalysis
            
            % (5) wrap to pi to bound the phase angle differences between 0 and pi:
            % delta_phi = wrapToPi(phi_ch1 - phi_ch2);
            euler_delta_phi = exp(1i*(phi_ch1 - phi_ch2));
        
            % mean vector (in complex space)
            mean_complex_vector = mean(euler_delta_phi); % magnitude of the mean
            ISPC = abs(mean_complex_vector); % it is equal to ISPC computed bellow
            
            % euler representation of angles
            delta_phi = angle(euler_delta_phi);
            PLI = abs(mean(sign(delta_phi)));% it is equal to PLI computed bellow
        
            figure()
            subplot(321)
            plot(1:numel(amplitude_ch1), [amplitude_ch1; amplitude_ch2])
            axis tight
            title('Instantaneous amplitude')
            % xlabel('Samples')
            ylabel('Voltage (\muV)')
            legend(eegChannelsName{electrode_pairs(ee,1)}, ...
                eegChannelsName{electrode_pairs(ee,2)})
            
            subplot(322)
            plot(1:numel(phi_ch1), amplitude_ch1-amplitude_ch2, 'Color', [0.5 0.5 0.5])
            axis tight
            title('Amplitude difference')
            % xlabel('Samples')
            ylabel('Voltage (\muV)')
            
            subplot(323)
            plot(1:numel(phi_ch1), [phi_ch1; phi_ch2])
            axis tight
            title('Instantaneous phase')
            xlabel('Samples')
            ylabel('Phase angle (radians)')
            legend(eegChannelsName{electrode_pairs(ee,1)}, ...
                eegChannelsName{electrode_pairs(ee,2)})
            
            subplot(324)
            plot(1:numel(phi_ch1), phi_ch1-phi_ch2, 'Color', [0.5 0.5 0.5])
            axis tight
            title('Phase difference')
            xlabel('Samples')
            ylabel('Phase angle (rad)')
            
            subplot(325)
            phase2plot_ch1 = [zeros(1,numel(phi_ch1))' phi_ch1']';
            polarplot(phase2plot_ch1(:), repmat([0 1],1,numel(phi_ch1))), hold on
            phase2plot_ch2 = [zeros(1,numel(phi_ch2))' phi_ch2']';
            polarplot(phase2plot_ch2(:), repmat([0 1],1,numel(phi_ch2))), hold off
            title('Phase angles from two channels')
            
            subplot(326)
            phase2plot = [zeros(1,numel(phi_ch1))' (phi_ch1 - phi_ch2)']';
            polarplot(phase2plot(:), repmat([0 1],1,numel(phi_ch1)), ...
                'Color', [0.5 0.5 0.5]), hold on
            polarplot([0 angle(mean_complex_vector)], [0 ISPC], ...
                'm', 'Linewidth', 5)
            title('Phase angle differences between two channels')
            
        end
        
        % length of mean vector corresponding to phase_synchronization or ISPC
        % (intersite phase clustering, Cohen2014, chapter 26)
        tic
        Sxy = z_ch1.*conj(z_ch2);
        Sxy_time = toc;
        tic
        ISPC = abs(mean(exp(1i*angle(Sxy))));
        
        count_matrices = count_matrices+1;
        count_measures = count_measures+1;
        undirect_measures_matrix(electrode_pairs(ee,1), electrode_pairs(ee,2), count_measures) = ISPC;
        undirect_measures_matrix(electrode_pairs(ee,2), electrode_pairs(ee,1), count_measures) = ISPC;
        undirect_measures_vector(ee, count_measures) = ISPC;
        undirect_measures_names(count_measures) = {['ISPC_' freq_band_names{ff}]};
        comp_time_vector(ee, count_measures) = toc+Sxy_time+get_analytic_time+get_elec_time+init_time;
        
        
        % (6) get Phase Lag Index (Niso2013):
        tic
        Sxyi = imag(Sxy);
        Sxyi_time = toc;
        tic
        PLI = abs(mean(sign(Sxyi)));
        
        count_matrices = count_matrices+1;
        count_measures = count_measures+1;
        undirect_measures_matrix(electrode_pairs(ee,1), electrode_pairs(ee,2), count_measures) = PLI;
        undirect_measures_matrix(electrode_pairs(ee,2), electrode_pairs(ee,1), count_measures) = PLI;
        undirect_measures_vector(ee, count_measures) = PLI;
        undirect_measures_names(count_measures) = {['PLI_' freq_band_names{ff}]};
        comp_time_vector(ee, count_measures) = toc+Sxyi_time+Sxy_time+get_analytic_time+get_elec_time+init_time;
        
        
        % weighted phase-lag index
        % WPLI is an extension of the PLI in which angle differences are
        % weighted according to their distance from the real axis
        tic
        WPLI = abs(mean(abs(Sxyi).*sign(Sxyi)))./mean(abs(Sxyi));
        
        count_matrices = count_matrices+1;
        count_measures = count_measures+1;
        undirect_measures_matrix(electrode_pairs(ee,1), electrode_pairs(ee,2), count_measures) = WPLI;
        undirect_measures_matrix(electrode_pairs(ee,2), electrode_pairs(ee,1), count_measures) = WPLI;
        undirect_measures_vector(ee, count_measures) = WPLI;
        undirect_measures_names(count_measures) = {['WPLI_' freq_band_names{ff}]};
        comp_time_vector(ee, count_measures) = toc+Sxyi_time+Sxy_time+get_analytic_time+get_elec_time+init_time;
        
        
        % debiased weighted phase-lag index (shortcut, as implemented in fieldtrip)
        tic
        imagsum      = sum(Sxyi);
        imagsumW     = sum(abs(Sxyi));
        debiasfactor = sum(Sxyi.^2);
        dWPLI  = (imagsum.^2 - debiasfactor)./(imagsumW.^2 - debiasfactor);
        
        count_matrices = count_matrices+1;
        count_measures = count_measures+1;
        undirect_measures_matrix(electrode_pairs(ee,1), electrode_pairs(ee,2), count_measures) = dWPLI;
        undirect_measures_matrix(electrode_pairs(ee,2), electrode_pairs(ee,1), count_measures) = dWPLI;
        undirect_measures_vector(ee, count_measures) = dWPLI;
        undirect_measures_names(count_measures) = {['dWPLI_' freq_band_names{ff}]};
        comp_time_vector(ee, count_measures) = toc+Sxyi_time+Sxy_time+get_analytic_time+get_elec_time+init_time;
        
        
        
        %% get Spearman's correlation coefficient for instantaneous power
        % correlations
        tic
        power_ch1 = amplitude_ch1.^2;
        power_ch2 = amplitude_ch2.^2;
        get_power_time = toc;
        tic
        RHO = corr(power_ch1', power_ch2', 'Type', 'Spearman');
        % covariance of two variables, scaled by the variance of each variable
        % assumption: data are normally distributed
        
        count_matrices = count_matrices+1;
        count_measures = count_measures+1;
        undirect_measures_matrix(electrode_pairs(ee,1), electrode_pairs(ee,2), count_measures) = RHO;
        undirect_measures_matrix(electrode_pairs(ee,2), electrode_pairs(ee,1), count_measures) = RHO;
        undirect_measures_vector(ee, count_measures) = RHO;
        undirect_measures_names(count_measures) = {['corr_coef_Spearman_power_' freq_band_names{ff}]};
        comp_time_vector(ee, count_measures) = toc+get_power_time+get_analytic_time+get_elec_time+init_time;
        
        
        %% get normalized cross-correlation (Mormann2003a, Mormann2005,
        % Mirowski2009)
        tic
        max_r = max(xcorr(power_ch1, power_ch2, maxlag, 'normalized'));
        
        count_matrices = count_matrices+1;
        count_measures = count_measures+1;
        undirect_measures_matrix(electrode_pairs(ee,1), electrode_pairs(ee,2), count_measures) = max_r;
        undirect_measures_matrix(electrode_pairs(ee,2), electrode_pairs(ee,1), count_measures) = max_r;
        undirect_measures_vector(ee, count_measures) = max_r;
        undirect_measures_names(count_measures) = {['corr_max_xcorr_power_' freq_band_names{ff}]};
        comp_time_vector(ee, count_measures) = toc+get_power_time+get_analytic_time+get_elec_time+init_time;
        
    end
    
    %% Directed Transfer Function
    
    % dtf = 1;
    % ind = 10;
    % undirect_measures_matrix(electrode_pairs(ee,1), electrode_pairs(ee,2), ind) = dtf;
    % undirect_measures_matrix(electrode_pairs(ee,2), electrode_pairs(ee,1), ind) = dtf;
    % undirect_measures_vector(ee, ind) = dtf;
    
    %% Partial Directed Coherence
    
    % pdc = 1;
    % ind = 11;
    % undirect_measures_matrix(electrode_pairs(ee,1), electrode_pairs(ee,2), ind) = pdc;
    % undirect_measures_matrix(electrode_pairs(ee,2), electrode_pairs(ee,1), ind) = pdc;
    % undirect_measures_vector(ee, ind) = pdc;
    
end

undirect_measures_names = undirect_measures_names(1:count_measures);
undirect_measures_vector = undirect_measures_vector(:,1:count_measures);
undirect_measures_matrix = undirect_measures_matrix(:,:,1:count_measures);
comp_time_vector = mean(comp_time_vector(:,1:count_measures));

%% INTERPRETATION OF MEASURES: ********************************************
% (1) PLI
% 0<=PLI<=1. (Kida2016, Niso2013, Cohen2014)
% "If all phase angle differences are on one side of the imaginary
% axis, the PLI will be high. If half of the phase angle differences
% are positive and half are negative with respect to the imaginary
% axis, the PLI will be zero." Cohen2014

% (2) ISPC = Phase Locking Value = Mean Phase Coherence (Niso2013, Cohen2014)
% 0<=ISPC<=1. Niso2013, Cohen2014
% values of phase synchronization above 0.5 is likely to be influenced
% by volume conduction (Cohen2014); this is a nondirectional measure of
% phase synchronization

% (3) 0<=WPLI<=1. 0 --> no synchronization, 1 --> synchronization (Niso2013)



%% Extract measures from the matrices computed above

% get circular omega complexity measure:
COC_features = zeros(n_bands,1);
COC_feat_names = cell(n_bands,1);
get_circ_corr = strfind(undirect_measures_names, 'Circular_correlation');
log_ind_circ_corr = ~cellfun(@isempty,get_circ_corr);
ind_circ_corr = find(log_ind_circ_corr);

for ff = 1:n_bands
    circ_corr_matrix = squeeze(undirect_measures_matrix(:,:,ind_circ_corr(ff)));
    eigen_values = eig(circ_corr_matrix);
    eigen_values = eigen_values/sum(eigen_values);
    om = sum(eigen_values.*log2(eigen_values));
    COC_features(ff) = abs(1 + om/log2(n_chans)); % COC measure value
    COC_feat_names(ff) = {['COC_' freq_band_names{ff}]};
end
% Baboukani2019: "COC varies between 0 and 1, where 0 indicates no phase
% synchrony between channels"

features = COC_features;
feature_names = COC_feat_names;

% threshold to exclude edges with non-significant interaction strengths
% threshold_vec = [0.8, 0.8, 0.8, 0.5, 0.5, 0.5, 0];
% VanMierlo2013: "We derived the 99th percentile for each connection
% between each pair of electrode contacts and used this as the threshold
% for that specific connection. In this way we expect to see only 1% of
% connections during the analysis of an iEEG epoch." even though they
% analysed effective connectivity

%%
plotFigureGraph = 0;


for mm = 1:numel(undirect_measures_names)
    % compute network measures for each connectivity matrix:
    % disp(undirect_measures_names{mm})
    
    % the Spearman's coefficient and the circular correlation matrix range
    % from -1 to 1, however to compute graph measures such as the
    % clustering coefficient, that requires all weights to be between 0 and
    % 1, I then took the absolute of the connectivity matrices
    % in terms of interpretation: for Spearman I'm observing both positive
    % and negative correlations, for the circular correlation I'm also
    % observing correlations that might be different in the phase (see
    % interpretation above)
    connectivity_matrix = abs(squeeze(undirect_measures_matrix(:,:,mm)));
    connectivity_vector = abs(undirect_measures_vector(:,mm));
    
    check_measure_PLI = strfind(undirect_measures_names{mm}, 'PLI');
    if ~isempty(check_measure_PLI) && check_measure_PLI==1 % PLI
        % as the PLI of each channel with itself is zero,
        % assign the diagonal with zeros:
        connectivity_matrix(1:n_chans+1:end) = 0;
        threshold = 0.5;
    elseif ~isempty(check_measure_PLI) && check_measure_PLI>1 % WPLI and dWPLI
        % as the WPLI and dWPLI of each channel with itself is NaN,
        % assign the diagonal with NaNs:
        connectivity_matrix(1:n_chans+1:end) = NaN;
        threshold = 0.5;
    end
    
    check_measure_corr = strfind(undirect_measures_names{mm}, 'corr');
    if ~isempty(check_measure_corr)
        if check_measure_corr(1)==1
            % for correlations that range between zero and one
            threshold = 0.8;
        else
            % for circular correlation
            threshold = 0.3;
        end
    end
    
    check_measure_ISPC = strfind(undirect_measures_names{mm}, 'ISPC');
    if ~isempty(check_measure_ISPC) && check_measure_ISPC==1
        threshold = 0.5;
    end
    
    connectivity_matrix = abs(connectivity_matrix);
    [MD, MS, CPL, GE, MCC, MBC, WGCC, T, M, A] = getNetworkUndirectedMeasures(...
        connectivity_matrix, connectivity_vector, eegChannelsName, ...
        electrode_pairs, threshold);
    
    features = [features; MD; MS; CPL; GE; MCC; MBC; WGCC; T; M; A];
    feature_names = [feature_names; strcat(undirect_measures_names{mm}, ...
        {'_MD'; '_MS'; '_CPL'; '_GE'; '_MCC'; '_MBC'; '_WGCC'; '_T'; ...
        '_M'; '_A'})];
    
    
    % Threshold defined by Abbas2021: 0.3 (with no given reason)
    if plotFigureGraph
        plotConnectivityMatrix(connectivity_matrix, connectivity_vector, ...
            eegChannelsName, electrode_pairs, threshold, ...
            undirect_measures_names{mm})
    end
    
end

% export_fig(gcf, fullfile(cd, 'FiguresPaper', 'connectivity_matrix_ISPC_pat_402_seiz1.pdf'), ...
%     '-painters', '-transparent')
    

%%

% get Phase Slope Index which is a directed measure:
% PSI measures whether the slope of the phase lag is consistently positive 
% or negative over several adjacent frequency bins (computed from the 
% Fourier transform). The sign of the slope indicates whether the net 
% connectivity flows from region A to B or the reverse. By specifying the 
% frequency bands we make the phase-slope index a frequency-band-specific 
% phase-based measure of directed connectivity. (Cohen2014)

th_PSI = 2;
PSI_features = [];
PSI_feat_names = [];
comp_time_PSI = [];
for ff = 1:n_bands
    tic
    fc1 = freq_bands(ff,1);
    fc2 = freq_bands(ff,2);
        
    name_connectivity = ['PSI_' freq_band_names{ff}];
    
    PSI_band = data2psiX(win_all_chan_data, fs, [fc1 fc2]);
    comp_time_PSI = [comp_time_PSI, toc];
    
    [MS, CPL, GE, MCIC, MCOC, MBC, M, A] = getNetworkDirectedMeasures( ...
        PSI_band, eegChannelsName, th_PSI);


    if plotFigureGraph
        plotConnectivityMatrix(PSI_band, [], eegChannelsName, ...
            electrode_pairs, th_PSI, name_connectivity)
    end
    
    PSI_features = [PSI_features; MS; CPL; GE; MCIC; MCOC; MBC; M; A];
    PSI_feat_names = [PSI_feat_names; strcat(name_connectivity, ...
        {'_MS'; '_CPL'; '_GE'; '_MCIC'; '_MCOC'; '_MBC'; '_M'; '_A'})];
end

comp_time_vector = [comp_time_vector, comp_time_PSI];


features = [features; PSI_features];
feature_names = [feature_names; PSI_feat_names];
comp_time_names = [undirect_measures_names; strcat('PSI_', freq_band_names)'];


% https://rcweb.dartmouth.edu/~mvdm/wiki/doku.php?id=analysis:course-w16:week11


end