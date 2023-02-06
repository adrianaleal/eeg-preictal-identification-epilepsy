function [energy_detail_levels, energy_approximation_levels] = ...
    discr_wavelet_trans(input_signal, mother_wav, dec_level, fs, ...
    disp_band_frequencies, plotFigure)


% Input:
% - input_signal
% - Mother wavelet: mother_wav
% - Decomposition levels: dec_level
% - Sampling frequency: fs
% - Flag do display the obtained frequency bands in command window: disp_band_frequencies


% cd = [];
% upsampled_cd = [];

% dec_level = 6;
% mother_wav = 'db4'

% Perform a multi-level decomposition of the signal using mother wavelet
[C,L] = wavedec(input_signal, dec_level, mother_wav);


if plotFigure
    length_upsampled_cds = zeros(dec_level,1);
    length_upsampled_cas = length_upsampled_cds;
end

energy_detail_levels = zeros(dec_level,1);
energy_approximation_levels = energy_detail_levels;

for n = 1:dec_level
    
    % wrcoef reconstructs the coefficients of a 1-D signal
    eval(['A{' num2str(n) '} = wrcoef(''a'',C,L,mother_wav,n);']) % 1-D detail coefficients
    eval(['D{' num2str(n) '} = wrcoef(''d'',C,L,mother_wav,n);']) % 1-D approximation coefficients
    
    
    
    % Extraction of the aproximation and details coeficients
    % eval(['ca{' int2str(n) '} = appcoef(C,L,mother_wav,n);']); % 1-D approximation coefficients
    % eval(['cd{' int2str(n) '} = detcoef(C,L,n);']); % 1-D detail coefficients
    
    % Upsampling approximation components:
    % eval(['upsampled_ca{' int2str(n) '} = upsample(ca{' int2str(n) '}, 2^n );']);
    
    % Upsampling detail components:
    % eval(['upsampled_cd{' int2str(n) '} = upsample(cd{' int2str(n) '}, 2^n );']);
    
    % CHECK THAT UPSAMPLED IS NOT CORRECT
    % figure()
    % subplot(211)
    % plot(A{1})
    % subplot(212)
    % plot(ca{1})
    % plot(upsampled_ca{1})
    
    % signal energy (^2):
    detail_coefficients = eval(['D{' int2str(n) '}']);
    detail_coefficients_squared = detail_coefficients.^2;
    energy_detail_levels(n) = sum(detail_coefficients_squared)/ ...
        length(detail_coefficients_squared);
    
    
    approximation_coefficients = eval(['A{' int2str(n) '}']);
    approximation_coefficients_squared = approximation_coefficients.^2;
    energy_approximation_levels(n) = sum(approximation_coefficients_squared)/ ...
        length(approximation_coefficients_squared);
    
end

%% display a table with the corresponding band frequencies in each level

if disp_band_frequencies==1
    fn_1st_filter = fs/2;

    cA_mat = zeros(dec_level,2);
    cD_mat = zeros(dec_level,2);
    cA_name = cell(dec_level,1);
    cD_name = cell(dec_level,1);
    for ii = 1:dec_level
        cA_mat(ii,:) = [0, round(fn_1st_filter/(2^ii)*100)/100];
        cA_name(ii) = {['cA' num2str(ii)]};
        cD_mat(ii,:) = [round(fn_1st_filter/(2^ii)*100)/100, round(fn_1st_filter/(2^(ii-1))*100)/100];
        cD_name(ii) = {['cD' num2str(ii)]};
    end

    band_frequencies = table(cA_name,cA_mat,cD_name,cD_mat)

end
%%

if plotFigure

    % maxL = max(length_upsampled_cds);
    % decomposition_coefficients = zeros(dec_level,maxL);
    

    % for n = 1:dec_level
    %     dif = maxL - (length(eval(['upsampled_cd{' int2str(n) '}'])));
    %     add = zeros(1,dif);
    %     eval(['upsampled_cd{' int2str(n) '} = horzcat(upsampled_cd{' int2str(n) '},add);']);
    %     decomposition_coefficients = eval(['upsampled_cd{' int2str(n) '}']);
    % end

    time = 1/fs:1/fs:length(D{1})/fs;
    
    figure()
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    subplot(621)
    plot(time,D{1},'Color',[.7 .2 .37]),
    legend('D1'), axis tight
    title('Discrete wavelet transform detail coefficients')
    subplot(623)
    plot(time,D{2},'Color',[.07 .5 .09])
    legend('D2'), axis tight
    subplot(625)
    plot(time,D{3},'m')
    legend('D3'), axis tight
    subplot(627)
    plot(time,D{4},'Color',[0.28 0.55 0.96])
    legend('D4'), axis tight
    subplot(629)
    plot(time,D{5},'Color',[1 .5 0])
    legend('D5'), axis tight
    xlabel('Time (s)')

    subplot(622)
    plot(time,A{1},'Color',[.7 .2 .37])
    title('Discrete wavelet transform approximation coefficients')
    legend('A1'), axis tight
    subplot(624)
    plot(time,A{2},'Color',[.07 .5 .09])
    legend('A2'), axis tight
    subplot(626)
    plot(time,A{3},'m'),hold on
    legend('A3'), axis tight
    subplot(628)
    plot(time,A{4},'Color',[0.28 0.55 0.96])
    legend('A4'), axis tight
    subplot(6,2,10)
    plot(time,A{5},'Color',[1 .5 0])
    legend('A5'), axis tight
    xlabel('Time (s)')%, ylabel('\textbf{Amplitude (mV)}')

end



end


