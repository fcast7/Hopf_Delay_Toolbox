function Hopf_delays_feature_extraction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Whole-brain Hopf Bifurcation Model with time delay:
%
%                   Spectral and dynamical features on simulated data
%
%  1) Input: Define Simulation Parameter (a, K, MD) 
%
%  2) Input:  Load simulated time series (Zsave);
%     Output: Spectral features of simulated time series in the Whole-Brain Hopf Model
%             to save as Model_Spectral_Features.mat (uncomment fo measuring also Spectral
%            Entropy and Lempel-Ziv Complexity)
%
%  Joana Cabral 2017 joana.cabral@psych.ox.ac.uk
%  Francesca Castaldo francesca.castaldo.20@ucl.ac.uk
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath(genpath('C:\Users\fcast\OneDrive - University College London\CLUSTER\PROJECT[AAL]\New_Simu'));
% 1) Define simulation Parameters

a=-5; %bifurcation parameter of your simulations 
N=90; %nodes parcellation
dt_save=2e-3;
f=40; %Hz intrinsic frequency 
MD=0:1:10; % Reange of Mean Delay in ms
MD=MD.*1e-3; % Mean Delay in seconds
fbins=1000;
freqZ = (0:fbins-1)/(dt_save*fbins);

expK=1.8:0.1:2;
% [rows, columns] = size(expK);
% expK = reshape(sort(expK(:), 'descend'), [columns, rows])';
K=10.^(expK);

Sync                 = zeros(length(K),length(MD));
Meta                 = zeros(length(K),length(MD));
PeakFMean            = zeros(length(K),length(MD));
PeakFGlobal          = zeros(length(K),length(MD));
% SpecEntropy_Global = zeros(length(K), length(MD));
% lzcomplexity       = zeros(length(K), length(MD));

% 2) Generate spectral features given simulated time series

% figure('color', 'w') %only if one wants to see the simulations in the whole parameter space

for g=1:length(K)
    for d=1:length(MD)
        k=K(g);
        md=MD(d);
        
        disp(['Now K=' num2str(k) ', mean Delay = ' num2str(md*1e3) 'ms'])
        
        K_label=num2str(log10(k));
        ind_p=find(K_label=='.');
        if numel(ind_p)
            K_label(ind_p)='p';
        end
        
        %Load your simulations (depending on the parameter you would like
        %to investigate)
%         cd('C:\Users\fcast\OneDrive - University College London\CLUSTER\PROJECT[AAL]\New_Simu')
%         load(['Hopf_Simu_K1E' K_label '_MD' num2str(MD*1e3)],'Zsave','dt_save')
        load(['a_Remote_K1E' K_label '_MD_' num2str(md*1e3) 'a' num2str(a)],'Zsave') %for cluster new simulations
        
        % Detect the peak Frequency in the Fourier Transform of all areas
        Fourier_Complex = fft(Zsave,fbins,2); %% Fourier of Z (complex) in 2nd dimension
       
        Fourier_Global=abs(mean(Fourier_Complex)).^2;
        [~, Imax]=max(Fourier_Global);
        PeakFGlobal(g,d)=freqZ(Imax);
        
        % Evaluate Order Stability around the Global peak frequency
        % Band-pass filter the signals Z around the peak frequency of the
        % ensemble
       
        for n=1:90
            Zb(n,:)=bandpasshopf(Zsave(n,:),[max(0.1,freqZ(Imax)-1) freqZ(Imax)+1],1/dt_save);
            Zb(n,:)=angle(hilbert(Zb(n,:)));
            
        end
        OP =abs(mean(exp(1i*(Zb)),1));
        
        Sync(g,d)=mean(OP); % measure of global synchronization (mean of OP)
        Meta(g,d)=std(OP); %  how much OP fluctuates in time (std of OP)
        
        % Probability distribution (total sum is 1)
        PSD_G=Fourier_Global/sum(Fourier_Global);
        
        
        % Average of Node PSD
        Fourier_Mean=mean(abs(Fourier_Complex).^2);
        [~, Imax]=max(Fourier_Mean);
        PeakFMean(g,d)=freqZ(Imax);
        
        % Probability distribution (total sum is 1)
        PSD_M=Fourier_Mean/sum(Fourier_Mean);
        
        %Spectral Entropy method 1
%         PSD=PSD_G/sum(PSD_G+1e-12);
%         logPSD = log2(PSD_G + 1e-12);
%         SpecEntropy_Global(g,d) = -sum(PSD_G.*logPSD)/log2(length(PSD_G)); %Normalized spectral entropy
        
        %         %Spectral Entropy Method 2 + lempel ziv complexity
        %         for n=1:90
        %           SE = pentropy(real(Zsave(n,:)), 1/dt_save, 'Instantaneous', false);
        %           LZC= LZComp_median(Zsave(n,:));
        %         end
        %         se(g,d) =SE;
        %         lzcomplexity(g,d) = LZC;
        %
        
        
        
        % Uncomment to generate a plot for simulations in all the parameter space
        % ( This is what fills the figure created at the beginning)
        
        %         subplot(length(g),length(d),d+(g-1)*numel(MD))
        %         plot(PSD_G)
        %         xlim ([0 50])
        %         axis off
        
        
        
    end
end
cd('C:\Users\fcast\OneDrive - University College London\CLUSTER\PROJECT[AAL]\New_Simu')
save([('Model_Spectral_Features')],'PeakFGlobal', 'PeakFMean','Sync', 'Meta') %,'lzcomplexity','SpecEntropy_Global')
