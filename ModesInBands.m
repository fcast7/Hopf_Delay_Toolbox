
addpath('Hopf_Delay_Toolbox');

% Simulations
simulation_files={'d4_HCP_Sim_Cluster_K1E1_MD0','d4_HCP_Sim_Cluster_K1E-1_MD3','Independent_AAL_HCP_Simu_K1E1_MD3','d4_HCP_Sim_Cluster_K1E1p7_MD3','AAL_HCP_Simu_K1E1_MD10'};
simulation_names={'No delays','Weak K','Intermediate K','Strong K','Long delays'};

N = 90;
Order=[1:2:N N:-2:2];

% Frequency bands to analyze (in Hz)
delta = [0.5 4];
theta = [4 8];
alpha = [8 13];
beta = [13 30];
gamma = [30 60];

n_simu=length(simulation_names);

figure
lim=0.002;

for simu=1:n_simu
    
    disp(simulation_names{simu})
    
    load(simulation_files{simu},'Zsave','dt_save')
    Zsave  = Zsave./(5*std(Zsave(:))); % Check if we remove this
    Zsave=Zsave(Order,:);
    Zdelta = zeros(size(Zsave));
    Ztheta = zeros(size(Zsave));
    Zalpha = zeros(size(Zsave));
    Zbeta  = zeros(size(Zsave));
    Zgamma  = zeros(size(Zsave));
    
    for n=1:90
        Zdelta(n,:) = bandpasshopf(Zsave(n,:), delta , 1/dt_save);
        Ztheta(n,:) = bandpasshopf(Zsave(n,:), theta , 1/dt_save);
        Zalpha(n,:) = bandpasshopf(Zsave(n,:), alpha , 1/dt_save);
        Zbeta(n,:)  = bandpasshopf(Zsave(n,:), beta  , 1/dt_save);
        Zgamma(n,:) = bandpasshopf(Zsave(n,:), gamma  , 1/dt_save);
    end
    
    % Remove the first and last seconds after band pass filtering
    Zdelta(:,[1:1/dt_save end-1/dt_save:end])=[];
    Ztheta(:,[1:1/dt_save end-1/dt_save:end])=[];
    Zalpha(:,[1:1/dt_save end-1/dt_save:end])=[];
    Zbeta(:,[1:1/dt_save end-1/dt_save:end])=[];
    Zgamma(:,[1:1/dt_save end-1/dt_save:end])=[];  
    
    % Calculate amplitude envelope
    Env_Delta=abs(hilbert(Zdelta'))';
    Env_Theta=abs(hilbert(Ztheta'))';
    Env_Alpha=abs(hilbert(Zalpha'))';
    Env_Beta =abs(hilbert(Zbeta'))';
    Env_Gamma =abs(hilbert(Zgamma'))';
     
    
    % Envelope Covariance and Eigenvalues (Modes)    
    FC_Delta=cov(Env_Delta');
    FC_Theta=cov(Env_Theta');
    FC_Alpha=cov(Env_Alpha');
    FC_Beta=cov(Env_Beta');    
    FC_Gamma=cov(Env_Gamma');  
    

% If you want to use Correlation instead of covariance 
%     FC_Delta=corrcoef(Env_Delta');
%     FC_Theta=corrcoef(Env_Theta');
%     FC_Alpha=corrcoef(Env_Alpha');
%     FC_Beta=corrcoef(Env_Beta');    
%     FC_Gamma=corrcoef(Env_Gamma');  

    EigVal_Delta=sort(eig(FC_Delta),'descend');
    EigVal_Theta=sort(eig(FC_Theta),'descend');
    EigVal_Alpha=sort(eig(FC_Alpha),'descend');
    EigVal_Beta=sort(eig(FC_Beta),'descend');
    EigVal_Gamma=sort(eig(FC_Gamma),'descend');
    
    
    Isubdiag = find(tril(ones(90),-1));

    colormap(jet)
    
    subplot(4,n_simu,simu)
    lim=4*std(abs(FC_Beta(Isubdiag)));
    %lim=mean(diag(FC_Beta));
    imagesc(FC_Beta,[-lim lim])
    axis square
    xticks([])
    yticks([])
    title(simulation_names(simu),'FontSize',12,'FontName','Helvetica','Interpreter','none')
    colorbar
    
    subplot(4,n_simu,simu+1*n_simu)
    lim=4*std(abs(FC_Alpha(Isubdiag)));
    %lim=mean(diag(FC_Alpha));
    imagesc(FC_Alpha,[-lim lim])
    axis square
    xticks([])
    yticks([])
    colorbar
    
    subplot(4,n_simu,simu+2*n_simu)
    lim=4*std(abs(FC_Theta(Isubdiag)));
    %lim=mean(diag(FC_Theta));
    imagesc(FC_Theta,[-lim lim])
    axis square
    xticks([])
    yticks([])
    colorbar
    
    subplot(4,n_simu,simu+3*n_simu)
    lim=5*std(FC_Delta(Isubdiag));
    %lim=mean(diag(FC_Delta));
    imagesc(FC_Delta,[-lim lim])
    axis square
    xticks([])
    yticks([])
    colorbar
    
    
    if simu==1
        EigVal_Delta_Thres= EigVal_Delta(1);
        EigVal_Theta_Thres= EigVal_Theta(1);
        EigVal_Alpha_Thres= EigVal_Alpha(1);
        EigVal_Beta_Thres = EigVal_Beta(1);
        EigVal_Gamma_Thres= EigVal_Gamma(1);
    else
        N_Modes_Delta(simu)=sum(EigVal_Delta>EigVal_Delta_Thres);
        N_Modes_Theta(simu)=sum(EigVal_Theta>EigVal_Theta_Thres);
        N_Modes_Alpha(simu)=sum(EigVal_Alpha>EigVal_Alpha_Thres);
        N_Modes_Beta(simu) =sum(EigVal_Beta>EigVal_Beta_Thres);   
        N_Modes_Gamma(simu)=sum(EigVal_Gamma>EigVal_Gamma_Thres);
    
        [Spatial_Modes_Delta{simu}, Val_Delta]=eigs(cov(Env_Delta'),N_Modes_Delta(simu));
        [Spatial_Modes_Theta{simu}, Val_Theta]=eigs(cov(Env_Theta'),N_Modes_Theta(simu));
        [Spatial_Modes_Alpha{simu}, Val_Alpha]=eigs(cov(Env_Alpha'),N_Modes_Alpha(simu));   
        [Spatial_Modes_Beta{simu}, Val_Beta]=eigs(cov(Env_Beta'),N_Modes_Beta(simu)); 
        [Spatial_Modes_Gamma{simu}, Val_Gamma]=eigs(cov(Env_Gamma'),N_Modes_Gamma(simu)); 

        
        Val_Modes_Delta{simu}=diag(Val_Delta);
        Val_Modes_Theta{simu}=diag(Val_Theta);
        Val_Modes_Alpha{simu}=diag(Val_Alpha);
        Val_Modes_Beta{simu}=diag(Val_Beta);
        Val_Modes_Gamma{simu}=diag(Val_Gamma);
    end    

end

figure

n_collumns=1;

% Strong K

for simu=1:n_simu
    subplot(1,n_simu,simu)
    hold on
    barh(1,N_Modes_Delta(simu), 'Facecolor', [1 .8 .5])
    barh(2,N_Modes_Theta(simu), 'Facecolor', [1 .7 .7])
    barh(3,N_Modes_Alpha(simu), 'Facecolor', [.7 .7 1])
    barh(4,N_Modes_Beta(simu), 'Facecolor', [.7 1 .7])
    %barh(5,N_Modes_Gamma(simu), 'Facecolor', [.7 1 .7])
    yticks([1 2 3 4])
    yticklabels({'\delta','\theta', '\alpha', '\beta'})
    ylabel('Frequency bands','FontSize',12,'FontName','Helvetica')
    xlabel('Eigenvalues > baseline','FontSize',12,'FontName','Helvetica')
    title(simulation_names(simu),'FontSize',12,'FontName','Helvetica','Interpreter','none')
    ylim([0 5])
    xlim([0 28])
end



% %% Render modes
% 
% simu=3;
% 
% figure
% 
% % delta modes
% colormap(jet)
% for mode_num=1:size(Spatial_Modes_Delta{simu},2)
%     disp(mode_num)
%     subplot_tight(4,25,3*25+mode_num,0.001)
%     V=Spatial_Modes_Delta{simu}(:,mode_num);
%     Render_links_in_cortex(V,[1 .7 0])
%     title(['\lambda= ' num2str(Val_Modes_Delta{simu}(mode_num),'%.2e')],'FontSize',7)
% end
% 
% % theta modes
% colormap(jet)
% for mode_num=1:size(Spatial_Modes_Theta{simu},2)
%     disp(mode_num)
%     subplot_tight(4,25,2*25+mode_num,0.001)
%     V=Spatial_Modes_Theta{simu}(:,mode_num);
%     Render_links_in_cortex(V,'r')
%     title(['\lambda= ' num2str(Val_Modes_Theta{simu}(mode_num),'%.2e')],'FontSize',7)
% end
% 
% % alpha modes
% colormap(jet)
% for mode_num=1:size(Spatial_Modes_Alpha{simu},2)
%     disp(mode_num)
%     subplot_tight(4,25,25+mode_num,0.001)
%     V=Spatial_Modes_Alpha{simu}(:,mode_num);
%     Render_links_in_cortex(V,'b')
%     title(['\lambda= ' num2str(Val_Modes_Alpha{simu}(mode_num),'%.2e')],'FontSize',7)
% end
% 
% % beta modes
% colormap(jet)
% for mode_num=1:size(Spatial_Modes_Beta{simu},2)
%     disp(mode_num)
%     subplot_tight(4,25,mode_num,0.001)
%     V=Spatial_Modes_Beta{simu}(:,mode_num);
%     Render_links_in_cortex(V,'g')
%     title(['\lambda= ' num2str(Val_Modes_Beta{simu}(mode_num),'%.2e')],'FontSize',7)
% end
% 
% simu=4;
% 
% figure
% 
% % delta modes
% colormap(jet)
% for mode_num=1:size(Spatial_Modes_Delta{simu},2)
%     disp(mode_num)
%     subplot_tight(4,25,3*25+mode_num,0.001)
%     V=Spatial_Modes_Delta{simu}(:,mode_num);
%     Render_links_in_cortex(V,[1 .7 0])
%     title(['\lambda= ' num2str(Val_Modes_Delta{simu}(mode_num),'%.2e')],'FontSize',7)
% end
% 
% % theta modes
% colormap(jet)
% for mode_num=1:size(Spatial_Modes_Theta{simu},2)
%     disp(mode_num)
%     subplot_tight(4,25,2*25+mode_num,0.001)
%     V=Spatial_Modes_Theta{simu}(:,mode_num);
%     Render_links_in_cortex(V,'r')
%     title(['\lambda= ' num2str(Val_Modes_Theta{simu}(mode_num),'%.2e')],'FontSize',7)
% end
% 
% % alpha modes
% colormap(jet)
% for mode_num=1:size(Spatial_Modes_Alpha{simu},2)
%     disp(mode_num)
%     subplot_tight(4,25,25+mode_num,0.001)
%     V=Spatial_Modes_Alpha{simu}(:,mode_num);
%     Render_links_in_cortex(V,'b')
%     title(['\lambda= ' num2str(Val_Modes_Alpha{simu}(mode_num),'%.2e')],'FontSize',7)
% end
% 
% % beta modes
% colormap(jet)
% for mode_num=1:size(Spatial_Modes_Beta{simu},2)
%     disp(mode_num)
%     subplot_tight(4,25,mode_num,0.001)
%     V=Spatial_Modes_Beta{simu}(:,mode_num);
%     Render_links_in_cortex(V,'g')
%     title(['\lambda= ' num2str(Val_Modes_Beta{simu}(mode_num),'%.2e')],'FontSize',7)
% end
