%% To plot all the source simulated psd (takes some time)
% 
function plot_simu_psd

cd('C:\Users\fcast\OneDrive - University College London\CLUSTER\PROJECT[AAL]\New_Simu\32AAL\a_neg5')
load ('MEG_Fitting.mat');
fbins=1000;
dt_save=2e-3;
freqZ = (0:fbins-1)/(dt_save*fbins);
ind150Hz=find(freqZ==150); %range of interest
% [rows, columns] = size(expK);
% expK = reshape(sort(expK(:), 'descend'), [columns, rows])'; %reverse to plot in the right order 
K=10.^(expK);

figure;
% for g=-length(K):1
g=1:length(K);
[rows, columns] = size(g);
k= reshape(sort(g(:), 'descend'), [columns, rows])'; %reverse to plot in the right order

for g=1:length(K)

    for d=1:length(MD)
     
        subplot(length(expK),length(MD),d+(g-1)*numel(MD))
        plot(freqZ(1:ind150Hz),squeeze(PSD_Simu_Global(k(g),d,:)),'k','LineWidth',0.8)
        xlim ([0 80])
        axis off
    end 

end
