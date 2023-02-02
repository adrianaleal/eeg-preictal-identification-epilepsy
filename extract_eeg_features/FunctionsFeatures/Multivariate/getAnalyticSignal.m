function [z, hilbert_amplitude, hilbert_phase] = getAnalyticSignal( ...
    data2filter, fs, fc1, fc2)

% Get the complex valued analytic signal or the time-frequency
% decomposition


time = linspace(1,numel(data2filter)/fs,numel(data2filter));

% Remove the mean of the data before filtering
data2filter = detrend(data2filter);



%% first option: filter-Hilbert

% 5th-order butterworth filter
norder = 5;
fn = fs/2;
wc1 = fc1/fn; % (center_freq-filter_frequency_spread)
wc2 = fc2/fn; % (center_freq+filter_frequency_spread)

% rp = 3; % at most 3 dB of ripple 
% rs = 40; % at least 40 dB attenuation in the stopbands
% n_but = buttord(fc1/fn, fc2/fn, rp, rs)

[butterB, butterA] = butter(norder, [wc1, wc2] ,'bandpass');
butter_filter = filtfilt(butterB, butterA, data2filter);

% analytic complex signal
z = hilbert(butter_filter);

% get the instantaneous phase angles:
hilbert_phase = angle(z); % in radians

% get the instantaneous amplitudes:
% hilbert_amplitude = real(z);
hilbert_amplitude = butter_filter;


% figure()
% subplot(211)
% % plot real part of the filtered signal
% plot(time, hilbert_amplitude, 'b'), hold on
% plot(time, butter_filter, 'r--')
% xlabel('Time (s)'), ylabel('Amplitude (\muV)');
% legend({'Butterworth'})
% 
% % now plot phases
% subplot(212)
% plot(time, hilbert_phase,'r')
% xlabel('Time (s)'), ylabel('Phase angles (rad.)');
% legend({'Butterworth'})
% 
% figure()
% plot(time, data2filter), hold on
% plot(time, hilbert_amplitude, 'r')


%% second option: Morlet wavelet convolution


end

