fclose all; clear; close all; clc;
set(groot, 'defaultAxesFontSize',12,'defaultTextFontSize',12)
set(groot, 'defaultAxesTickLabelInterpreter','tex');
set(groot, 'defaultLegendInterpreter','tex');
set(groot, 'defaultTextInterpreter','tex')

% Required toolboxes:
% - Signal Processing Toolbox
% - Statistics and Machine Learning Toolbox
% - Wavelet Toolbox
% - Predictive Maintenance Toolbox

% add folder containing the functions to compute EEG features
addpath(genpath(fullfile(cd, 'FunctionsFeatures')))

% add utils folder
addpath(genpath(fullfile(cd, 'utils')))


n_chans = 19;
sampling_frequency = 256;

% load example of 5 second EEG window for all channels
load('example_5sec_eeg_window.mat')

% load channels' names
load('eeg_channels_name_example_5sec_eeg_window.mat')

% define the groups of features that are to be extracted
group_features2extract = {'univariate_linear'};
% options: 'univariate_linear', 'univariate_nonlinear', 'multivariate'


if any(strcmp(group_features2extract, 'univariate_linear') | ...
        strcmp(group_features2extract, 'univariate_nonlinear'))
    
    for eegChannelIndex = 1:n_chans
        
        get_5sec_win_data_chan = get_5sec_win_all_chan_data(eegChannelIndex,:);
        
        %% Linear univariate features
        if any(strcmp(group_features2extract, 'univariate_linear'))
            plotFigure = 1;
            [features_univariate_linear, features_names_univariate_linear, ...
                comp_times_univariate_linear, comp_times_names_univariate_linear] ...
                = univariate_linear_features(get_5sec_win_data_chan, ...
                sampling_frequency, plotFigure);
        end
        
        %% Nonlinear univariate features
        if any(strcmp(group_features2extract, 'univariate_nonlinear'))
            plotFigure = 1;
            [features_univariate_nonlinear, features_names_univariate_nonlinear, ...
                comp_times_univariate_nonlinear, comp_times_names_univariate_nonlinear] = ...
                univariate_nonlinear_features(get_5sec_win_data_chan, ...
                sampling_frequency, plotFigure);
        end
        
    end
end

if any(strcmp(group_features2extract, 'multivariate'))
   [features_multivariate, features_names_multivariate, ...
        comp_times_multivariate, comp_times_names_multivariate] = ...
        multivariate_features(get_5sec_win_all_chan_data, ...
        sampling_frequency, eeg_channels_name); 
end
