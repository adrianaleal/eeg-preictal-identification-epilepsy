function [features, feature_names, comp_times, comp_times_names] = ...
    univariate_linear_features(input_signal, fs, plotFigure)


% Input:
% - input_signal: signal to analyse
% - fs: sampling frequency
% - plotFigure: flag to plot figure

% Output:
% - features: univariate linear features
% - feature_names: names of features
% - comp_times: computation time for each group of features
% - comp_times_names: name of each group of features



%% Statistical moments (BouAssi2017, Mormann2007, Valderrama2012)
tic
std_vector = std(abs(input_signal));
kurtosis_vector = kurtosis(input_signal);
skewness_vector = skewness(input_signal);

mean_intensity_normalized = mean(abs(input_signal))/max(abs(input_signal));
mean_intensity = mean(abs(input_signal));
ct_statistical_moments = toc;

%% computation of the power spectrum

tic

% Hamming window
% win_length = 2^6;
[psd_signal, frequency_vector] = pwelch(input_signal, [], [], [], fs);
% [pxx,f] = pwelch(x, window, noverlap, nfft, fs)
% window default: Hamming window is used to obtain eight segments of x 
% with noverlap overlapping samples. = length(input_signal)/8 = 160
% noverlap default: 50% overlap between segments
% nfft default: max(256,2^nextpow2(length(window)))
% nfft = max(256,2^nextpow2(160)) = 256
% frequency resolution = 128 Hz / 256 FFT bins  = 0.5 Hz
% welch_overlap_samples = floor(n_samples * welch_overlap / 100);
% "The frequency resolution of a Fourier transform is defined by the 
% number of points in the time series. Therefore, the larger the time 
% segment, the more frequencies can be extracted, and thus the greater 
% the frequency resolution. The temporal resolution is always defined 
% by the data sampling rate and thus does not depend on the length of 
% the time segment." Cohen2014

ct_psd = toc;

%% Spectral band power features (BouAssi2017, Mormann2007, Valderrama2012)

tic
% psd in the frequency band waves
delta_band = [0.5, 4];
theta_band = [4, 8];
alpha_band = [8, 13];
beta_band = [13, 30];
gamma1_band = [30, 47];
gamma2_band = [53, 75];
gamma3_band = [75, 97];
gamma4_band = [103, 128];
% this frequency splitting in the gamma band should exclude powerline 
% interference from the data (BouAssi2017)

% calculating the power (area below curves)
total_power = trapz(abs(psd_signal));

% The relative spectral power is computed as the EEG signals contain much
% more power in the low-frequency range BouAssi2017

freq_band_names = {'delta'; 'theta'; 'alpha'; 'beta'; ...
    'gamma1'; 'gamma2'; 'gamma3'; 'gamma4'};

for ff = 1:numel(freq_band_names)
    
    
    eval([freq_band_names{ff} '_psd = psd_signal(intersect(find(frequency_vector>=' ...
        freq_band_names{ff} '_band(1)), find(frequency_vector<' freq_band_names{ff} '_band(2))));'])
    
    if total_power~=0 
        % compute the power (area below curves)
        eval([freq_band_names{ff} '_power = trapz(abs(' ...
            freq_band_names{ff} '_psd));'])
        
        % compute the relative spectral power
        eval([freq_band_names{ff} '_relative_power = ' freq_band_names{ff} ...
            '_power/total_power;'])
    else
        % compute the power (area below curves)
        eval([freq_band_names{ff} '_power = 0;'])
        
        % compute the relative spectral power
        eval([freq_band_names{ff} '_relative_power = 0;'])
    end
end


features_aux = [delta_power; theta_power; alpha_power; beta_power; ...
    gamma1_power; gamma2_power; gamma3_power];

% alpha peak
alpha_band_indexes = intersect(find(frequency_vector>=alpha_band(1)), ...
    find(frequency_vector<alpha_band(2)));
alpha_peak_idx = find(alpha_psd==max(alpha_psd));
alpha_peak_freq = frequency_vector(alpha_band_indexes(alpha_peak_idx));
alpha_peak_freq = alpha_peak_freq(1);


% ratio features
combs = nchoosek(1:numel(features_aux),2);
ratio_features_names = join(freq_band_names(combs), '_');

n_comb = length(combs);
ratio_bands_vector = zeros(n_comb,1);
n_ratio_bands = numel(features_aux);
count = 0;
for f=1:n_ratio_bands-1
    for g=f+1:n_ratio_bands
        count = count+1;
        ratio_bands_vector(count,:) = features_aux(f,:)./features_aux(g,:);
    end
end


beta2alpha_theta_ratio = beta_power/(alpha_power+theta_power);
theta2alpha_beta_ratio = theta_power/(alpha_power+beta_power);

% the mean frequency of a power spectral density (PSD) estimate, pxx. The 
% frequencies, f, correspond to the estimates in pxx:
mean_frequency = meanfreq(abs(psd_signal),frequency_vector);

ct_frequency_feat = ct_psd+toc;

%% Spectral edge power and frequency 
% (BouAssi2017, Mormann2007, Valderrama2012)
tic

[spectral_edge_frequency, spectral_edge_cum_power, spectral_edge_power] ...
    = spectral_edge(frequency_vector, psd_signal);

ct_SEF = ct_psd+toc; 

%% Hjorth parameters (BouAssi2017, Mormann2007, Valderrama2012)

tic
% Activity
activity = var(input_signal);
% Mobility
mobility = std(diff(input_signal))/std(input_signal);
% Complexity
complexity = std(diff(diff(input_signal)))./std(diff(input_signal))./mobility;

ct_hjorth = toc;

%% Decorrelation time (BouAssi2017, Valderrama2012, Direito2017, Teixeira2014m)

tic

autocorrelation_vector = xcorr(input_signal, 'unbiased');
% decorrelation time: the first zero crossing of the autocorrelation function
decorrelation_time = find(autocorrelation_vector(length(input_signal):end)<=0,1);

ct_dec_time = toc;

%% Accumulated energy or Long-term energy (BouAssi2017, Mormann2007)
% according to "Accumulated energy revisited" study there are no
% differences over interictal to motivate the used of this feature

%% Discrete wavelet transform (BouAssi2017, Faust2015)

disp('*******************************************************************')
disp('*************** DISCRETE WAVELET TRANSFORM ************************')
disp('*******************************************************************')
tic

mother_wav = 'db4';
dec_level = 5;
disp_band_frequencies = 1;
[energy_detail_levels, energy_approximation_levels] = discr_wavelet_trans( ...
    input_signal, mother_wav, dec_level, fs, disp_band_frequencies, ...
    plotFigure);

features_dwt = [energy_detail_levels; energy_approximation_levels(5)];
features_dwt_names = {'energy_D1'; 'energy_D2'; 'energy_D3'; ...
    'energy_D4'; 'energy_D5'; 'energy_A5'};

ct_wavelets = toc;

%%

features = [delta_power; theta_power; alpha_power; beta_power; ...
    gamma1_power; gamma2_power; gamma3_power; gamma4_power; ...
    delta_relative_power; theta_relative_power; alpha_relative_power; ...
    beta_relative_power; gamma1_relative_power; gamma2_relative_power; ...
    gamma3_relative_power; gamma4_relative_power; total_power; ...
    alpha_peak_freq; mean_frequency; ratio_bands_vector; ...
    beta2alpha_theta_ratio; theta2alpha_beta_ratio; ...
    mean_intensity_normalized; mean_intensity; ...
    std_vector; kurtosis_vector; skewness_vector; activity; mobility; ...
    complexity; spectral_edge_frequency; spectral_edge_cum_power; ...
    decorrelation_time; features_dwt];

feature_names = [{'Delta_power'; 'Theta_power'; 'Alpha_power'; ...
    'Beta_power'; 'Gamma1_power'; 'Gamma2_power'; 'Gamma3_power'; ...
    'Gamma4_power'; 'Relative_delta_power'; 'Relative_theta_power'; ...
    'Relative_alpha_power'; 'Relative_beta_power'; ...
    'Relative_gamma1_power'; 'Relative_gamma2_power'; ...
    'Relative_gamma3_power'; 'Relative_gamma4_power'; 'Total_power'; ...
    'Alpha_peak_frequency'; 'Mean_frequency'}; ...
    strcat({'Ratio_'}, ratio_features_names); ...
    {'Ratio_beta_over_alpha_theta'; 'Ratio_theta_over_alpha_beta'; ...
    'Normalized_mean_intensity'; 'Mean_intensity'; 'Std'; ...
    'Kurtosis'; 'Skewness'; 'Activity'; 'Mobility'; 'Complexity'; ...
    'Spectral_edge_frequency'; 'Spectral_edge_power'; ...
    'Decorrelation_time'}; features_dwt_names];


comp_times = [ct_statistical_moments; ct_frequency_feat; ct_SEF; ...
    ct_hjorth; ct_dec_time; ct_wavelets];


comp_times_names = {'ct_statistical_moments'; 'ct_frequency_feat'; ...
    'ct_SEF50'; 'ct_hjorth'; 'ct_dec_time'; 'ct_wavelets'};

%% plot figure

if plotFigure
    figure()
    set(gcf,'units','normalized','outerposition',[0 0 0.5 0.8])
    subplot(211)
    time = 1/fs:1/fs:length(input_signal)/fs;
    plot(time, input_signal)
    axis tight
    ylabel('uV')
    xlabel({'Time (s)', '(a)'})
    
    subplot(212)
    h1 = plot(frequency_vector, psd_signal, 'k');
    axis tight
    hold on
    ylabel('Power Spectral Density (uV^2/Hz)');
    % xlabel('Frequency (Hz)');
    
    legend_entries = {'Welch''s PSD estimate'};
    
    yrange = ylim;
    color_vec = [1 1 0; 1 0 0; 0 0 1; 0 0.5 0; 0.850 0.325 0.098; ...
        0.850 0.325 0.098; 0.850 0.325 0.098; 0.850 0.325 0.098];
    face_vec = [0.5 0.5 0.5 0.5 0.9 0.7 0.5 0.2];
    
    % plot the frequency bands:
    save_h = h1;
    
    for bb = 1:numel(freq_band_names)
        eval(['band_h = fill([' freq_band_names{bb} '_band(1) ' ...
            freq_band_names{bb} '_band(1) ' freq_band_names{bb} ...
            '_band(2) ' freq_band_names{bb} '_band(2)], [yrange(1) yrange(2) yrange(2) yrange(1)], color_vec(bb,:));'])
        band_h.FaceAlpha = face_vec(bb);
        band_h.EdgeAlpha = face_vec(bb);
        legend_entries = [legend_entries, ...
            {[upper(freq_band_names{bb}(1)) freq_band_names{bb}(2:end)]}];
        save_h = [save_h, band_h];
    end
    
    plot(frequency_vector, psd_signal, 'k')
    
    % plot the spectral edge frequency:
    h2 = plot(spectral_edge_frequency, spectral_edge_power, 'rv', ...
        'MarkerSize', 8, 'MarkerFaceColor', 'red');
    legend_entries = [legend_entries, {'f_{50}'}];
    
    % plot the alpha peak:
    h3 = plot(alpha_peak_freq, max(alpha_psd), 'bv', ...
        'MarkerSize', 8, 'MarkerFaceColor', 'blue');
    legend_entries = [legend_entries, {'Alpha peak frequency'}];
    
    
    axis tight
    xlabel({'Frequency (Hz)', '(b)'})
    
    hold off
    legend([save_h, h2, h3],legend_entries)
%     saveas(gcf, 'frequency_analysis.emf')
%     export_fig(gcf, '-dpdf', fullfile(cd, 'frequency_analysis.pdf'))
    
end

end