function Network_Parameter_Space

%  1) Input: load DEMO simu: Model_Spectral_Features.mat and create plots
%
%             
%     Output: Parameter space plots (Global Peak Frequency, Mean Peak Frequency, Metastability,
%            Spectral Entropy )
%
%  2) Theoretical Prediction by Niebur Paper 1991


% 1) Plot GlobalPeakF, MeanPeakF, Meta, Synch

% LOAD DEMO DATA: PeakFGlobal, PeakFMean, Sync, Meta parameter range of interest


load('Model_Spectral_Features.mat');

dt_save=2e-3;
f=40; %Hz
MD=0:1:10; % Reange of Mean Delay in ms
MD=MD.*1e-3; % Mean Delay in seconds
fbins=1000;
freqZ=(0:fbins-1)/(dt_save*fbins);

expK=-1:0.1:1.7;
% [rows, columns] = size(expK); %this is needed because of how K is measured in the simulations 
% expK = reshape(sort(expK(:), 'descend'), [columns, rows])';
K=10.^(expK); 

% Create figure
fig=figure ('color', 'w');
colormap(jet);
imagesc(MD*1e3,log10(K),PeakFGlobal)
title('Global Peak Frequency (Hz)','FontSize',20,'FontName','Helvetica')
axis xy
ylabel('Global Coupling K','FontSize',14,'FontName','Helvetica');
xlabel('Mean Conduction Delay (ms)','FontSize',14,'FontName','Helvetica')
yticklabels({'0.1','','1','','10',''});
colorbar
caxis([0 45])
h = colorbar;
ylabel(h, 'Frequency (Hz)','FontSize',10,'FontName','Helvetica')
saveas(fig,'PeakF.png')

fig=figure ('color', 'w');
colormap(jet)
imagesc(MD*1e3,log10(K),PeakFMean)
title('Mean Peak Frequency (Hz)','FontSize',20,'FontName','Helvetica')
axis xy
ylabel('Global Coupling K','FontSize',14,'FontName','Helvetica');
xlabel('Mean Conduction Delay (ms)','FontSize',14,'FontName','Helvetica')
yticklabels({'0.1','','1','','10',''});
colorbar
caxis([0 45])
h = colorbar;
ylabel(h, 'Frequency (Hz)','FontSize',10,'FontName','Helvetica')
saveas(fig,'MeanF.png')

fig=figure ('color', 'w');
colormap(jet)
imagesc(MD*1e3,log10(K),Sync)
axis xy
ylabel('Global Coupling K','FontSize',14,'FontName','Helvetica');
xlabel('Mean Conduction Delay (ms)','FontSize',14,'FontName','Helvetica')
yticklabels({'0.1','','1','','10',''});
colorbar
caxis([0 1])
title('Synchrony mean(KOP)','FontSize',20,'FontName','Helvetica');
saveas(fig,'Synch.png')

%
fig=figure ('color', 'w');
colormap(jet)
imagesc(MD*1e3,log10(K),Meta)
title('Metastability std(KOP)','FontSize',20,'FontName','Helvetica');
ylabel('Global Coupling K','FontSize',14,'FontName','Helvetica');
xlabel('Mean Conduction Delay (ms)','FontSize',14,'FontName','Helvetica')
yticklabels({'0.1','','1','','10',''});
axis xy
colorbar
saveas(fig,'Meta.png')
% caxis([0 0.2])
% saveas(fig,'Meta.png')
% subplot(2,3,3)
% imagesc(K,MD*1e3,PeakF')
% axis xy
% ylabel('Mean Delay (ms)')
% xlabel('Coupling K')
% colorbar
% title('Oscillation Frequency')
% colormap(jet)

%% If you also had measured Spectral Entropy in Function 1 
%  in DEMO simulations there is no Spectral Entropy calculation

% fig=figure ('color', 'w');
% colormap(jet)
% imagesc(MD*1e3,log10(K),SpecEntropy_Global)
% title('Spectral Entropy','FontSize',20,'FontName','Helvetica');
% ylabel('Global Coupling K','FontSize',14,'FontName','Helvetica');
% xlabel('Mean Delay (ms)','FontSize',14,'FontName','Helvetica')
% yticklabels({'0.1','','1','','10',''});
% axis xy
% colorbar
% caxis([0 1])
% % saveas(fig,'Spec_Entropy.png')

% subplot(2,3,3)
% imagesc(K,MD*1e3,PeakF')
% axis xy
% ylabel('Mean Delay (ms)')
% xlabel('Coupling K')
% colorbar
% title('Oscillation Frequency')
% colormap(jet)



