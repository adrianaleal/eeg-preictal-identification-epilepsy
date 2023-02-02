function [spectral_edge_frequency, spectral_edge_cum_power, ...
    spectral_edge_power] = spectral_edge(frequency_vector, psd_signal)

% Spectral edge power and frequency 
% (BouAssi2017, Mormann2007, Valderrama2012)
% lowest frequency up to which half of the total power is contained in a 
% signal. it is equivalent to the median frequency of a signal(Ouyang2020).
% In the context of seizure prediction, however, an x value of 50% has been
% successfully employed to describe the minimum frequency up to which 50% 
% of the overall power of the 0–40 Hz band is contained. The area under 
% spectral edge frequency is called spectral edge power;


% Input:
% - frequency_vector
% - psd_signal

% Output:
% - spectral_edge_frequency
% - cumulative sum of power at spectral edge frequency used as feature: 
%   spectral_edge_cum_power
% - power at spectral edge frequency used for ploting reasons: 
%   spectral_edge_cum_power


% get the indexes of the 0–40 Hz frequency band:
indexes_band = find(frequency_vector<=40);

% get the 0–40 Hz frequency and power vectors:
frequency_band = frequency_vector(indexes_band);
power_band = psd_signal(indexes_band);

% get the total power in that band:
total_spectral_power_band = sum(power_band);

% define the percentage of overall power:
power_threshold = 50;

% get the the corresponding value of power:
spectral_power_threshold = total_spectral_power_band*(power_threshold/100);

% get the cumulative power:
cumulative_power_band = cumsum(power_band);

% get the the corresponding value of frequency and power:
ind_spectral_edge = find(cumulative_power_band>=spectral_power_threshold,1);
spectral_edge_frequency = frequency_band(ind_spectral_edge);
spectral_edge_cum_power = cumulative_power_band(ind_spectral_edge);
spectral_edge_power = power_band(ind_spectral_edge);

end