clear all;close all; clc;
addpath(genpath('Codigo_nao_linear'))
addpath(genpath('hurst estimators'))

person={'P1','P2','P3'};

%%%%%%%%%%%%%% Por Regioes %%%%%%%%%%%%%%%
system={'FPZ','F3','FZ','F4','T7','C3','CZ','C4','T8','P3','PZ','P4','OZ'};
a=[1:13];

chan_name=system;

dificuldades={'dificuldade1','dificuldade2','dificuldade3'};

evento={'cruz1','texto','cruz2','codigo','cruz3'};

n_channels=length(a);

comp_time=[];

for p=1
    tic
    clear name
    mkdir(strcat('nonlinear_features\',person{p}))
    
    load(strcat('C:\Users\julio\OneDrive\Área de Trabalho\dados_eeglab\remocao_inicio_fim\inter_referenced\cont_remo\ICA\clean\divididos\regioes\',person{p},'.mat'))

    for k=1:13
        tic
        k
        mkdir(strcat('nonlinear_features\',person{p},'\',chan_name{a(k)}))
        for m=1:3
            mkdir(strcat('nonlinear_features\',person{p},'\',chan_name{a(k)},'\',dificuldades{m}))
            for n=1:4
                
                signal=double(name.(dificuldades{m}).(evento{n})(a(k),:));
                time=[1:1:length(signal)];
                fs=1000;
                
                % window size
                window_time=10;             
                overlap=0.8;
                non_overlap=0.2;
                
                %time vector with overlap
                time1=linspace(time(1),time(end)-fs*window_time,round((time(end)-fs*window_time*overlap)/(fs*window_time*non_overlap)));
                time2=linspace(time(window_time*1000),time(end),round((time(end)-fs*window_time*overlap)/(fs*window_time*non_overlap)));
                
                fractal_dimension_vector=linspace(1,0,length(time1));
                hurst_vector=linspace(1,0,length(time1));
                multifractal_vector=ones(11,length(time1));
                ApEn_vector=linspace(1,0,length(time1));
                lyapExpMat_vector=linspace(1,0,length(time1));
                corrDimMat_vector=linspace(1,0,length(time1));
                simplicity_vector=linspace(1,0,length(time1));
                
                %                     plotFigure=0;
                
                for i=1:length(time1)
                                        tic
                    beginning_index=uint64(time1(i));
                    ending_index=uint64(time2(i));
                    
                    signal_segment=signal(beginning_index:ending_index);
                    
                    %                     tic
                    %                                             disp('>> Compute Katz Fractal Dimension')
                    %                                             fractaldim=Katz_FD(signal_segment);
                    fractaldim=Higuchi_FD(signal_segment,100);
                    fractal_dimension_vector(i)=fractaldim;
                    %                     end
                    %                     figure;plot(fractal_dimension_vector)
                    
                    %                     end
                    %                     toc
                    
                    %                     tic
                    %                     disp('>> Compute Hurst Exponent')
                    hurst=hurst_exponent(signal_segment);
                    hurst_vector(i)=hurst;
                    %                     toc
                    
                    %                     tic
                    %                     disp('>> Compute Multi Fractal Dimension')
                    [dhpre,hpre] = dwtleader(signal_segment);
                    dhpre=flip(dhpre);
                    hpre=flip(hpre);
                    hmin=hpre(1);
                    dmin=dhpre(1);
                    hmax=hpre(11);
                    dmax=dhpre(11);
                    hpic=hpre(6);
                    dH=hmax-hmin;
                    dD=dmin-dmax;
                    dmin_dpic=dhpre(6)-dmin;
                    dpic_dmax=dhpre(6)-dmin;
                    hmin_hpic=hpic-hmin;
                    hpic_hmax=hmax-hpic;
                    
                    multifractal_vector(:,i)=[hmin;dmin;hmax;dmax;hpic;dH;dD;dmin_dpic;dpic_dmax;hmin_hpic;hpic_hmax];
                    %                     toc
                    
                    %% normalize signal
                    
                    sig2analyse = signal_segment; % signal_segment
                    
                    sig_norm = (sig2analyse-mean(sig2analyse))/std(sig2analyse); % VER LIVRO Varsavsky2011
                    
                    %% Embedding dimension, calculated with False Nearest Neighbors(FNN)
                    %                     tic
                    %                     disp('>> Phase Space Reconstruction')
                    [attractorMat, tauMat, eDimMat] = phaseSpaceReconstruction(sig_norm);
                    %                     toc
                    
                    %                     tic
                    %                     disp('>> Computing approximate entropy')
                    ApEn = approximateEntropy(sig_norm', 'Dimension', eDimMat, 'Lag',tauMat);
                    ApEn_vector(i)=ApEn;
                    %                    toc
                    
                    %                    tic
                    %                     disp('>> Compute Correlation Largest Lyapunov Exponent ')
                    lyapExpMat = lyapunovExponent(sig_norm, fs, 'Dimension', eDimMat, 'Lag',tauMat, 'ExpansionRange', [1 30]);
                    lyapExpMat_vector(i)= lyapExpMat;
                    %                     toc
                    
                    %                     tic
                    %                     disp('>> Compute Simplicity measure ')
                    emb_matrix=(attractorMat/sqrt(length(attractorMat)))'*(attractorMat/sqrt(length(attractorMat)));
                    eigvalues=sort(eig(emb_matrix),'descend');
                    eigvalues=eigvalues./sum(eigvalues);
                    H_simpl=-sum(eigvalues.*log10(eigvalues));
                    simp_value=1/(2^H_simpl);
                    simplicity_vector(i)=simp_value;
                    %                     toc
                    
                    %                     tic
                    %                      disp('>> Computing Correlation Dimension')
                    [corrDimMat,rRange,corInt] = correlationDimension(sig_norm,'Dimension', eDimMat, 'Lag',tauMat,'MinRadius',0.08);
                    corrDimMat_vector(i) = corrDimMat;
                                        a=toc
                    
                end
                
                
                names=["fractal_dimension_vector";...
                    "hurst_vector";...
                    "MF_minH";"MF_minD";"MF_maxH";"MF_maxD";"MF_Hpico";...
                    "MF_deltaH";"MF_deltaD";"MF_larg1";"MF_larg2";"MF_larg2por1";
                    "ApEn_vector";...
                    "lyapExpMat_vector";...
                    "corrDimMat_vector";...
                    "simplicity"];
                
                features=[fractal_dimension_vector;hurst_vector; ...
                    multifractal_vector;ApEn_vector;lyapExpMat_vector;...
                    corrDimMat_vector; simplicity_vector];
                
                resume={features,names};
                
            end
            
        end
        fprintf('\n')
        fprintf('[%.3f sec] >> Elapsed time for non linear features of one channel \n', toc);
        fprintf('\n')
    end
end

