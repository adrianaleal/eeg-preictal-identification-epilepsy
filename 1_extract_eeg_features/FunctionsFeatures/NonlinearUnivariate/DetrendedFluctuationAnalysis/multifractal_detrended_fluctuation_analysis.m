function [Ht,Htbin,Ph,Dh, ms_peak_abcissa, ms_width, ms_symmetry] = ...
    multifractal_detrended_fluctuation_analysis(input_signal, scale_range, ...
    convert2RandomWalk, plotFigure)

% DFA Detrended fluctuation analysis
% Method for determining the statistical self-affinity of a signal.
% Useful in revealing the extent of long-range correlations in time series
% Sources:
% - https://www.physionet.org/physiotools/dfa/
% - Ihlen2012 "Introduction to Multifractal Detrended Fluctuation Analysis
%   in Matlab" --> STRONGLY RECOMMEND READING
% - Bryce2012 "Revisiting detrended fluctuation analysis"
% - MoralesMartinez2021 "Detrended Fluctuation Analysis (MFDFA) approach
%   for multifractal analysis of precipitation"
% - Tang2020 "Seizure Prediction Using Multi-View Features and Improved
%   Convolutional Gated Recurrent Network"

%   Calculates the DFA of a signal and it's scaling exponent alpha.
%   Input:
%       - time: time (or x values of signal)
%       - input_signal: signal data (or y values of signal)
%       - scale_range: Range of segment size values to use for
%         calculating the alpha scaling exponent. The input should be a
%         matrix of the number of alpha ranges (rows) by two columns
%         corresponding to the minimum and maximum number of samples in a
%         segment, respectively.
%   Output:
%       - segment_size_vec: vector containing the different segments
%         (x-axis of DFA)
%       - n_segment_vec: vector containing the number of segments for each
%         segment size
%       - fn: DFA value for each segment size
%       - alpha_vec: Exponential scaling factors corresponding to the
%         ranges defined in scale_range matrix


% Limitations: at least 8000 data points should be used; monofractal method;
% (Bryce2012, Ihlen2012, MoralesMartinez2021)





%% CONSIDERATIONS ON ECG DATA:
% DFA is used to quantify the fractal scaling properties of heart rate
% signals of short interval. Describes the fractal-like correlation
% properties of R–R interval data (Hoshi2013)

% DFA plot is not strictly linear but rather consisted of two distinct
% regions of different slopes separated at a break point suggesting that
% there is a short range scaling exponent (alpha1) over periods of 4 to 11
% beats (or 4 to 13), and a long-range exponent (alpah2), over longer periods
% (larger than 11 beats) (Hoshi2013)

%%
N = length(input_signal);

if nargin == 1
    scale_range = [10 N/10];
    % Ihlen2012: "A maximum segment size below 1/10 of the sample size of
    % the time series will provide at least 10 segments in the computation
    % of Fq.
    
    plotFigure = 0;
end

if nargin == 2
    convert2RandomWalk = 1;
end

if nargin == 3
    plotFigure = 0;
end

%%



% define the order of the polynomial fit
m = 1; % if the polinomial trend is linear
% m = 2; % if the polinomial trend is quadratic
% m = 3; % if the polinomial trend is cubic




% get the minimum number of samples in a segment
min_segment_size = min(scale_range(:));
% Ihlen2012: "The minimum segment size larger than 10 samples is a "rule of
% tumb" for the computation of RMS."


% get the maximum number of samples in a segment
max_segment_size = max(scale_range(:));
% Ihlen2012: "A maximum segment size below 1/10 of the sample size of the
% time series will provide at least 10 segments in the computation of F."

% define the number of divisions:
n_scales = 20;

exponents = linspace(log2(min_segment_size),log2(max_segment_size),n_scales);
% exponents = linspace(log2(min_segment_size),log2(max_segment_size),n_scales);
segment_size_vec = unique(round(2.^exponents));
n_scales = length(segment_size_vec); % update the number of divisions
% Ihlen2012: "It is favorable to have a equal spacing between scales when
% they are represented in the log-log DFA."


% Convert a noise like time series into a random walk like time series:
if convert2RandomWalk==1
    rw_input_signal = cumsum(input_signal - mean(input_signal));
else
    rw_input_signal = input_signal;
end
% rw_input_signal = input_signal;
% OR integrate the signal without mean
% OR given a bounded time series of length N, integration or summation
% first converts this into an unbounded process
% OR rw_input_signal is called cumulative sum or profile



scale_small = 7:2:17;

%% get F(q==0)

% analyse multifractality:
q = linspace(-5,5,101);

if plotFigure
    Fq = NaN(numel(q), n_scales);
end

Fq0 = NaN(1, n_scales);

for n = 1:n_scales
    segment_size = segment_size_vec(n); % size of the segments
    n_segment = floor(N/segment_size);
    RMS = zeros(n_segment,1);
    
    for jj = 1:n_segment
        ind_start = ((jj-1)*segment_size)+1;
        ind_stop = jj*segment_size;
        indexes_segment = ind_start:ind_stop;
        signal_segment = rw_input_signal(indexes_segment);
        
        C = polyfit(indexes_segment, signal_segment, m);
        segment_trend = polyval(C,indexes_segment);
        residual_variation = signal_segment-segment_trend;
        RMS(jj) = sqrt(mean(residual_variation.^2));
    end
    
    if plotFigure
        % Compute the overall fluctuation or overall RMS OR scaling function
        % for each q-value
        for nq = 1:numel(q)
            qRMS = RMS.^q(nq);
            % overall q-order RMS:
            Fq(nq,n) = mean(qRMS).^(1/q(nq));
        end
        % because 1/0 goes to infinit:
        Fq(q==0,n) = exp(0.5*mean(log(RMS.^2)));
    end
    
    % because 1/0 goes to infinit:
    Fq0(n) = exp(0.5*mean(log(RMS.^2)));
end

Cq0 = polyfit(log10(segment_size_vec),log10(Fq0),m);
Regfitq0 = polyval(Cq0,log10(scale_small));
Hq0 = Cq0(1);

%% DFA

halfmax = floor(max(scale_small)/2);
Time_index = halfmax+1:N-halfmax;
Ht = zeros(numel(scale_small),numel(Time_index));

for n = 1:numel(scale_small)
    
    halfseg = floor(scale_small(n)/2);
    RMS = zeros(N,1);
    
    for jj = Time_index
        
        index = (jj-halfseg:jj+halfseg)';
        signal_segment = (rw_input_signal(index))';
        
        % (1) Fit a polynomial trend to each segment
        if m==1
            % OPTION1: lowest computational effort for a linear fit
            x = [ones(length(index), 1), index];
            slope = x\signal_segment; % x/y
            yn = x * slope;
            segment_trend = yn;
        else
            % OPTION2: for fitting polynomials of higher order
            % get the polinomial coefficients C used to create the
            % polynomial trend segments_trend(:, jj):
            C = polyfit(index, signal_segment, m);
            segment_trend = polyval(C, index);
            % OR:
            % segments_trend(:, jj) = indexes_segments(:, jj).*C(1)+C(2);
        end
        
        % (2) Compute the local fluctuation or RMS around a trend
        % segments_trend(:, jj), for each segment
        residual_variation = signal_segment-segment_trend;
        RMS(jj) = sqrt(mean(residual_variation.^2));
        
    end
    
    RMSt = RMS(Time_index);
    resRMS = Regfitq0(n)-log10(RMSt);
    logscale = log10(length(Time_index))-log10(scale_small(n));
    Ht(n,:) = resRMS./logscale+Hq0;
    
end


Ht_row = Ht(:);
BinNumb = round(sqrt(length(Ht_row)));
[freq,Htbin] = hist(Ht_row,BinNumb);
Ph = freq./sum(freq);
Ph_norm = Ph./max(Ph);
Dh = 1-(log(Ph_norm)./log(mean(diff(Htbin))));

% FEATURE ESTRACTED FROM THE MULTIFRACTAL SPECTRUM (Tang2020):
% Tang2020: "The abscissa value of the apex, which is the most
% prominent index in the spectrum to characterize the temporal
% structure of each channel, is used to form the time-domain feature
% matrix of the EEG segment."
ms_peak_abcissa = mean(Htbin(Ph==max(Ph)));
ms_width = Htbin(end)-Htbin(1);

% fitting a parabola to Dh ()
C = polyfit(Htbin(~isinf(Dh)), Dh(~isinf(Dh)), 2);
parabola = polyval(C, Htbin);
ms_symmetry = C(2);

% Sikdar2018: "For symmetric spectrum, C(2) is zero and for asymmetric 
% shape, its value can be positive for left skewed and negative for right 
% skewed spectrum."
if plotFigure
    
    
    % This plot is intended to verify the existence of a multifractal
    % behaviour in the input_signal under analysis (according to Ihlen2012)
    % Namely, if the figure shows that the slopes of the regression lines
    % are q-dependent we are in the presence of a multifractal time series.
    % The difference between the q-order RMS for positive and negative q's
    % are more visual apparent at the small segment sizes compared to
    % the large segment sizes
    
    % Both the monofractal and whitenoise time series have no periods with
    % small and large fluctuation (the regression lines are parallels)
    
    colors_vec = lines(7);
    
    figure()
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    subplot(2,6,1:3)
    Hq = zeros(1,numel(q));
    q_plot = [-5,-3,-1,0,1,3,5];
    save_leg_handle = zeros(1, numel(q_plot));
    indexes_q_plot = ismember(q,q_plot);
    count = 0;
    for nq = 1:numel(q)
        
        
        Cq = polyfit(log10(segment_size_vec),log10(Fq(nq,:)),m);
        Regfit = polyval(Cq, log10(segment_size_vec));
        
        % Get alpha line
        alpha_line = 10.^(Regfit);
        
        if ismember(q(nq),q_plot)
            count = count+1;
            loglog(segment_size_vec, Fq(nq,:), 'ko', 'MarkerSize', 5, ...
                'MarkerFaceColor', colors_vec(count,:), ...
                'Color', colors_vec(count,:)); hold on
       
            lh = loglog(segment_size_vec, alpha_line, ...
                'LineWidth', 2, 'Color', colors_vec(count,:));
            save_leg_handle(count) = lh;
        end
        
        Hq(nq) = Cq(1);
        
    end
    grid on
    axis tight
    hold off
    xlabel('log10(n)'); ylabel('log10[F(n)]');
    set(gca,'XTick', 2.^(1:15)); % Set ticks at powers of two
    % if the
    title('Multifractal time series')
    legend(save_leg_handle, strcat('q = ', cellstr(num2str(q_plot'))), ...
        'Location', 'southeast')
    
    subplot(2,6,4:6)
    plot(q_plot, Hq(indexes_q_plot), 'k'), hold on
    scatter(q_plot, Hq(indexes_q_plot), 60, colors_vec, 'filled') 
    hold off
    xlabel('q'), ylabel('q-order Hurst exponent')
    grid on, axis tight
    
    subplot(2,6,7:8)
    tq = Hq.*q-1;
    hq = diff(tq)./(q(2)-q(1));
    Dq = (q(1:end-1).*hq)-tq(1:end-1);    
    plot(hq, Dq, 'k')
    xlabel('hq'), ylabel('Dq')
    grid on, axis tight
    title('Multifractal spectrum of Dq and hq')
    
    subplot(2,6,9:10)
    plot(Htbin, Ph, 'k'), hold on
    plot(ms_peak_abcissa.*ones(2,1), [min(Ph),max(Ph)], ...
        'Color', [0.749 0 0.749]), hold off
    xlabel('Ht'), ylabel('Ph')
    grid on, axis tight
    title('The probability distribution Ph of the local Hurst exponents Ht')
    
    subplot(2,6, 11:12)
    plot(Htbin, Dh, '*-'), hold on
    plot(Htbin, parabola)
    plot(ms_peak_abcissa.*ones(2,1), [min(Dh(~isinf(Dh))),max(Dh)], ...
        'Color', [0.749 0 0.749]), hold off
    grid on, axis tight
    xlabel('Ht'), ylabel('Dh')
    title('Multifractal spectrum Dh estimated from distribution Ph')
    
end


end