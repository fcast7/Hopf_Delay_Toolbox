% function MOM_analysis

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to detect MOMs in different frequency bands for a selected point
% in the parameter space
%  
%  1. Load DEMO simulations (the one used to generate plots in the paper)
%     and Define threshold power from a baseline scenario.
%     Resonance is detected when power increases 5 STD above baseline
%  2. For a selected working point (SIMU), define the time points when the power
%     is above threshold in bands Delta, Theta, Alpha and Beta
%  3.  Measure MOMs size, Duration and Occupancy. Plot.
%  4.  Figure signals overtime (for each parcellated area) 

%  Needs: bandpasshopf.m, conver_back_to_time.m, subplot_tight.m
%
% July 2021
% Joana Cabral and Francesca Castaldo
% joanacabral@med.uminho.pt
% francesca.castaldo.20@ucl.ac.uk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a=-5;
% Frequency bands to analyze (in Hz)
delta = [0.5 4];
theta = [4 8];
alpha = [8 13];
beta = [13 30];
low_pass=30;

% 1- SIMULATION FILE TO READ and DEFINE BASELINE
% addpath('Hopf_Delay_Toolbox');
addpath(genpath('C:\Users\fcast\OneDrive - University College London\CLUSTER\PROJECT[AAL]\New_Simu'));
simulation_file={'a_Remote_K1E1_MD_0a-5.mat','a_Remote_K1E-1_MD_2a-5.mat','a_Remote_K1E1p1_MD_2a-5.mat','a_Remote_K1E1p7_MD_6a-5.mat','a_Remote_K1E1_MD_10a-5.mat'};
% simulation_file={'d4_HCP_Sim_Cluster_K1E1_MD0','d4_HCP_Sim_Cluster_K1E-1_MD3','Independent_AAL_HCP_Simu_K1E1_MD3','d4_HCP_Sim_Cluster_K1E1p7_MD3','AAL_HCP_Simu_K1E1_MD10'};
simulation_names={'No delays','Weak K','Intermediate K','Strong K','Long delays'};

% BASELINE:  Load data in 0 delay case, intermediate coupling to define the baseline
% load('d4_HCP_Sim_Cluster_K1E1_MD0.mat'); %a=-5;
load('a_Remote_K1E1_MD_0a-5.mat')

N = size(Zsave,1);
Order=[1:2:N N:-2:2];
Zsave  = Zsave./(5*std(Zsave(:)));
Zsave=Zsave(Order,:);
Zdelta = zeros(size(Zsave));
Ztheta = zeros(size(Zsave));
Zalpha = zeros(size(Zsave));
Zbeta  = zeros(size(Zsave));
Zfilt  = zeros(size(Zsave));

for n=1:90
    Zdelta(n,:) = bandpasshopf(Zsave(n,:), delta , 1/dt_save);
    Ztheta(n,:) = bandpasshopf(Zsave(n,:), theta , 1/dt_save);
    Zalpha(n,:) = bandpasshopf(Zsave(n,:), alpha , 1/dt_save);
    Zbeta(n,:)  = bandpasshopf(Zsave(n,:), beta  , 1/dt_save);
    Zfilt(n,:)=   bandpasshopf(Zsave(n,:),[0.01 low_pass],1/dt_save);
end
clear Zsave

% Remove the first and last seconds after band pass filtering
Zdelta(:,[1:1/dt_save end-1/dt_save:end])=[];
Ztheta(:,[1:1/dt_save end-1/dt_save:end])=[];
Zalpha(:,[1:1/dt_save end-1/dt_save:end])=[];
Zbeta(:,[1:1/dt_save end-1/dt_save:end])=[];
Zfilt(:,[1:1/dt_save end-1/dt_save:end])=[];

% Define thresholds as 5*STD of the power in each band in each area
deltathr=5*std(Zdelta,[],2);
thetathr=5*std(Ztheta,[],2);
alphathr=5*std(Zalpha,[],2);
betathr=5*std(Zbeta,[],2);
filtthr=5*std(Zfilt,[],2);
clear Zalpha Zbeta Ztheta Zfilt


%% 2 - Load data in the working point and filter into bands
% Change the "simu" to look at other scenario outside the optimal working
% point 

simu=2;
load(simulation_file{simu})
N = size(Zsave,1);
Order=[1:2:N N:-2:2];

Zsave=Zsave./(5*std(Zsave(:)));
Zsave=Zsave(Order,:);
Zdelta=zeros(size(Zsave));
Ztheta=zeros(size(Zsave));
Zalpha=zeros(size(Zsave));
Zbeta=zeros(size(Zsave));
Zfilt  = zeros(size(Zsave));

for n=1:90
    Zdelta(n,:) = bandpasshopf(Zsave(n,:),delta,1/dt_save);
    Ztheta(n,:) = bandpasshopf(Zsave(n,:),theta,1/dt_save);
    Zalpha(n,:) = bandpasshopf(Zsave(n,:),alpha,1/dt_save);
    Zbeta(n,:)  = bandpasshopf(Zsave(n,:),beta,1/dt_save);
    Zfilt(n,:)  = bandpasshopf(Zsave(n,:),[0.01 low_pass],1/dt_save);
end
clear Zsave

% Remove the first and last seconds after band passing
Zdelta(:,[1:1/dt_save end-1/dt_save:end])=[];
Ztheta(:,[1:1/dt_save end-1/dt_save:end])=[];
Zalpha(:,[1:1/dt_save end-1/dt_save:end])=[];
Zbeta(:,[1:1/dt_save end-1/dt_save:end])=[];
Zfilt(:,[1:1/dt_save end-1/dt_save:end])=[];

% Calculate amplitude envelope
Env_Delta =abs(hilbert(Zdelta'))';
Env_Theta =abs(hilbert(Ztheta'))';
Env_Alpha =abs(hilbert(Zalpha'))';
Env_Beta  =abs(hilbert(Zbeta'))';
Env_Zfilt =abs(hilbert(Zfilt'))';

clear Zdelta Zalpha Zbeta Ztheta

% Detect formation of a MOMs in each band
T_Delta =zeros(size(Env_Delta));
T_Theta =zeros(size(Env_Theta));
T_Alpha =zeros(size(Env_Alpha));
T_Beta  =zeros(size(Env_Beta));

for n=1:N
    T_Delta(n,:) =Env_Delta(n,:) >deltathr(n);
    T_Theta(n,:) =Env_Theta(n,:) >thetathr(n);
    T_Alpha(n,:) =Env_Alpha(n,:) >alphathr(n);
    T_Beta(n,:)  =Env_Beta(n,:)  >betathr(n);
end


clear Env_Delta Env_Theta Env_Alpha Env_Beta


%% 3 - Calculate MOM durations and size, Coalition size and Occupancy

% Delta MOM Durations
MOM_Durations=[];
for n=1:N
    % Detect switches in and out of this state
    a=find(diff(T_Delta(n,:))==1); %on
    b=find(diff(T_Delta(n,:))==-1); %off
    
    % We discard the cases where state sarts or ends ON
    if length(b)>length(a)
        b(1)=[];
    elseif length(a)>length(b)
        a(end)=[];
    elseif  ~isempty(a) && ~isempty(b) && a(1)>b(1)
        b(1)=[];
        a(end)=[];
    end
    if ~isempty(a) && ~isempty(b)
        MOM_Durations=[MOM_Durations b-a];
    end
end
Delta_Mean_Duration= mean(MOM_Durations)*dt_save;
Delta_std_Duration= std(MOM_Durations)*dt_save;


% Theta MOM Durations
MOM_Durations=[];
for n=1:N
    % Detect switches in and out of this state
    a=find(diff(T_Theta(n,:))==1); %on
    b=find(diff(T_Theta(n,:))==-1); %off
    
    % We discard the cases where state sarts or ends ON
    if length(b)>length(a)
        b(1)=[];
    elseif length(a)>length(b)
        a(end)=[];
    elseif  ~isempty(a) && ~isempty(b) && a(1)>b(1)
        b(1)=[];
        a(end)=[];
    end
    if ~isempty(a) && ~isempty(b)
        MOM_Durations=[MOM_Durations b-a];
    end
end
Theta_Mean_Duration= mean(MOM_Durations)*dt_save;
Theta_std_Duration= std(MOM_Durations)*dt_save;

% Alpha MOM Durations
MOM_Durations=[];
for n=1:N
    a=find(diff(T_Alpha(n,:))==1); %on
    b=find(diff(T_Alpha(n,:))==-1); %off
    
    if length(b)>length(a)
        b(1)=[];
    elseif length(a)>length(b)
        a(end)=[];
    elseif  ~isempty(a) && ~isempty(b) && a(1)>b(1)
        b(1)=[];
        a(end)=[];
    end
    if ~isempty(a) && ~isempty(b)
        MOM_Durations=[MOM_Durations b-a];
    end
end
Alpha_Mean_Duration= mean(MOM_Durations)*dt_save;
Alpha_std_Duration= std(MOM_Durations)*dt_save;

% Beta MOMs
MOM_Durations=[];
for n=1:N
    a=find(diff(T_Beta(n,:))==1); %on
    b=find(diff(T_Beta(n,:))==-1); %off
    if length(b)>length(a)
        b(1)=[];
    elseif length(a)>length(b)
        a(end)=[];
    elseif  ~isempty(a) && ~isempty(b) && a(1)>b(1)
        b(1)=[];
        a(end)=[];
    end
    if ~isempty(a) && ~isempty(b)
        MOM_Durations=[MOM_Durations b-a];
    end
end
Beta_Mean_Duration= mean(MOM_Durations)*dt_save;
Beta_std_Duration= std(MOM_Durations)*dt_save;

% 4 - Coalition size (N of elements contributing to the emergence of MOMs)

Delta_coalition_members=sum(T_Delta);
Delta_coalition_members=Delta_coalition_members(Delta_coalition_members>0);
Delta_coalition_members_mean=mean(Delta_coalition_members);
Delta_coalition_members_std=std(Delta_coalition_members);
clear Delta_coalition_members

Theta_coalition_members=sum(T_Theta);
Theta_coalition_members=Theta_coalition_members(Theta_coalition_members>0);
Theta_coalition_members_mean=mean(Theta_coalition_members);
Theta_coalition_members_std=std(Theta_coalition_members);
clear Theta_coalition_members

Alpha_coalition_members=sum(T_Alpha);
Alpha_coalition_members=Alpha_coalition_members(Alpha_coalition_members>0);
Alpha_coalition_members_mean=mean(Alpha_coalition_members);
Alpha_coalition_members_std=std(Alpha_coalition_members);

Beta_coalition_members=sum(T_Beta);
Beta_coalition_members=Beta_coalition_members(Beta_coalition_members>0);
Beta_coalition_members_mean=mean(Beta_coalition_members);
Beta_coalition_members_std=std(Beta_coalition_members);

clear Alpha_coalition_members Beta_coalition_members

% MOM Occupancy

Delta_Occ=sum(T_Delta(:))/numel(T_Delta);
Theta_Occ=sum(T_Theta(:))/numel(T_Theta);
Alpha_Occ=sum(T_Alpha(:))/numel(T_Alpha);
Beta_Occ = sum(T_Beta(:))/numel(T_Beta);

% Plot: FIGURE MEAN DURATION, SIZE, OCCUPANCY

% barplot
figure
subplot(1,3,1,'Linewidth',1,'Fontsize',14)
hold on
bar(1,Delta_Mean_Duration, 'Facecolor', [1 .8 .5],'Linewidth',1)
bar(2,Theta_Mean_Duration, 'Facecolor', [1 .7 .7],'Linewidth',1)
bar(3,Alpha_Mean_Duration, 'Facecolor', [.7 .7 1],'Linewidth',1)
bar(4,Beta_Mean_Duration, 'Facecolor', [.7 1 .7],'Linewidth',1)
ylim([0 1.2])
hold on
errorbar([Delta_Mean_Duration Theta_Mean_Duration Alpha_Mean_Duration Beta_Mean_Duration],[Delta_Mean_Duration Theta_std_Duration Alpha_std_Duration Beta_std_Duration], 'LineStyle', 'none')
xticks([1 2 3 4])
xticklabels({'\delta','\theta', '\alpha', '\beta'})
%xlabel('Frequency Bands','FontSize',14,'FontName','Helvetica')
ylabel('Duration (sec)','FontSize',14,'FontName','Helvetica')
title(simulation_file{simu},'FontSize',14,'FontName','Helvetica','Interpreter','none')

subplot(1,3,2,'Linewidth',1,'Fontsize',14)
hold on
bar(1,Delta_coalition_members_mean, 'Facecolor', [1 .8 .5],'Linewidth',1)
bar(2,Theta_coalition_members_mean, 'Facecolor', [1 .7 .7],'Linewidth',1)
bar(3,Alpha_coalition_members_mean, 'Facecolor', [.7 .7 1],'Linewidth',1)
bar(4,Beta_coalition_members_mean, 'Facecolor', [.7 1 .7],'Linewidth',1)
hold on
errorbar([Delta_coalition_members_mean Theta_coalition_members_mean Alpha_coalition_members_mean Beta_coalition_members_mean],[Delta_coalition_members_std Theta_coalition_members_std Alpha_coalition_members_std Beta_coalition_members_std], 'LineStyle', 'none')
xticks([1 2 3 4])
xticklabels({'\delta','\theta', '\alpha', '\beta'})
%xlabel('Frequency Bands','FontSize',14,'FontName','Helvetica')
ylabel('MOM size','FontSize',14,'FontName','Helvetica')
%title('Working point Delay=3ms, K=10','FontSize',18,'FontName','Helvetica')
ylim([0 70])

subplot(1,3,3,'Linewidth',1,'Fontsize',14)
hold on
bar(1,Delta_Occ, 'Facecolor', [1 .8 .5],'Linewidth',1)
bar(2,Theta_Occ, 'Facecolor', [1 .7 .7],'Linewidth',1)
bar(3,Alpha_Occ, 'Facecolor', [.7 .7 1],'Linewidth',1)
bar(4,Beta_Occ, 'Facecolor', [.7 1 .7],'Linewidth',1)
hold on
xticks([1 2 3 4])
xticklabels({'\delta','\theta', '\alpha', '\beta'})
%xlabel('Frequency Bands','FontSize',14,'FontName','Helvetica')
ylabel('MOM Occupancy','FontSize',14,'FontName','Helvetica')
%title('Working point Delay=3ms, K=10','FontSize',18,'FontName','Helvetica')
ylim([0 1])

%% 4- Figure Signals over time

Time_to_plot=25;
%Zfilt=Zfilt(:,1:Time_to_plot/dt_save);

Phase_filt=angle(hilbert(Zfilt'))';
OP=abs(mean(exp(1i*Phase_filt)));

load AAL_labels.mat label90

figure('Color','white')
hold on

% Delta patches
for n=1:N
    u=0;
    y=[];    x=[];
    for tx=find(T_Delta(n,1:Time_to_plot/dt_save))
        u=u+1;
        y(:,u)= [n-0.5 n-0.5  n+.5 n+.5];
        x(:,u) = [tx-1 tx tx tx-1];
    end
    p=patch(x.*dt_save,y,'y');
    set(p,'LineStyle','none','FaceColor','y','FaceAlpha',0.3);
end
ylim([0 N+1])

% Theta patches
for n=1:N
    u=0;
    y=[];    x=[];
    for tx=find(T_Theta(n,1:Time_to_plot/dt_save))
        u=u+1;
        y(:,u)= [n-0.5 n-0.5  n+.5 n+.5];
        x(:,u) = [tx-1 tx tx tx-1];
    end
    p=patch(x.*dt_save,y,'r');
    set(p,'LineStyle','none','FaceColor','r','FaceAlpha',0.3);
end
ylim([0 N+1])

% Alpha patches
for n=1:N
    u=0;
    y=[];    x=[];
    for tx=find(T_Alpha(n,1:Time_to_plot/dt_save))
        u=u+1;
        y(:,u)= [n-0.5 n-0.5  n+.5 n+.5];
        x(:,u) = [tx-1 tx tx tx-1];
    end
    p=patch(x.*dt_save,y,'b');
    set(p,'LineStyle','none','FaceColor','b','FaceAlpha',0.3);
end

for n=1:N
    u=0;
    y=[];
    x=[];
    for tx=find(T_Beta(n,1:Time_to_plot/dt_save))
        u=u+1;
        y(:,u)= [n-0.5 n-0.5  n+.5 n+.5];
        x(:,u) = [tx-1 tx tx tx-1];
    end
    p=patch(x.*dt_save,y,'g');
    set(p,'LineStyle','none','FaceColor','g','FaceAlpha',0.3);
end

plot(0:dt_save:(length(Zfilt)-1)*dt_save,(1:N)'.*ones(size(Zfilt))+(Zfilt./(2*filtthr)),'k')
title(['Simulated signal in the 90 units  coupled with K=' num2str(K) ' and MD= ' num2str(MD*1000) ' ms, filtered below ' num2str(low_pass) 'Hz'],'FontSize',14,'FontName','Helvetica')
xlabel('Time (seconds)','FontSize',18,'FontName','Helvetica')
ylabel('Units representing anatomically-defined cortical and subcortical brain areas','FontSize',18,'FontName','Helvetica')
% xlim([0 Interval])
box off
set(gca,'YTick',1:N,'Fontsize',8)
set(gca,'YTickLabel',[])
%set(gca,'YTickLabel',label90(Order,:))
ylim([0 N+1])
xlim([0 Time_to_plot])

[cc p]=corrcoef(sum(Env_Zfilt),OP)

OP=OP(1:Time_to_plot/dt_save);
Env_Zfilt=Env_Zfilt(:,1:Time_to_plot/dt_save);

figure
[AX,H1,H2]=plotyy(0:dt_save:(length(OP)-1)*dt_save,mean(Env_Zfilt),0:dt_save:(length(OP)-1)*dt_save,OP);
set(AX(2),'ytick',0:0.5:1,'Fontsize',16)
set(AX(1),'ytick',[0.1 0.2],'ylim',[0.05 0.25],'Fontsize',16)
xlabel('Time (seconds)','Fontsize',16)
legend({'Mean Amplitude Envelope','Kuramoto Order Parameter'},'Orientation','horizontal','Fontsize',16)
box off

%% Video MOMs over time

T_timepoints=1:10:5000;

Brain_Mask=niftiread('MNI152_T1_2mm_brain_mask.nii');
scortex=smooth3(Brain_Mask);
clear Brain_Mask

% Normalize Zfilt between 0 and 1 for the renderings

Zfilt=Zfilt/(2*std(Zfilt(:)));
Zfilt(Zfilt>1)=1;
Zfilt(Zfilt<-1)=-1;
Zfilt=(Zfilt+1)/2;

load AAL_cog_MNI2mm Parcels_cog
MNI_coord=Parcels_cog;
clear aal_cog

% PLOT SPHERES IN THE LOCATION OF EACH AREA
a=3;
[x,y,z] = sphere;
x=a*x;
y=a*y;
z=a*z;

figure('color','w')

% videoMOMs = VideoWriter(['VideoMOMs_' simulation_file '.mp4'],'MPEG-4'); %create the video object
% videoMOMs.FrameRate = 25;
% videoMOMs.Quality = 25;
% open(videoMOMs); %open the file for writing


hold on
% First plot a transparent cortex
cortexpatch=patch(isosurface(scortex,0.1), 'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'none','FaceAlpha', 0.1);
reducepatch(cortexpatch,0.1);
isonormals(scortex,cortexpatch);

% Then plot one sphere in each area
for n=1:N
    s(n)=surf(x+MNI_coord(n,2), y+MNI_coord(n,1),z+MNI_coord(n,3),'FaceColor',[0.9 .9 .9],'EdgeColor','none','FaceAlpha',0.2);
end

axis equal
material dull
lighting gouraud
daspect([1 1 1])
h = camlight('headlight');
xlim([5 105])
ylim([5 85])
zlim([0 80])
axis off

view(0,20)
%view(-90,90) % Top view

for t=1:length(T_timepoints)
    Delta_color=T_Delta(:,T_timepoints(t));
    Alpha_color=T_Alpha(:,T_timepoints(t));
    Beta_color=T_Beta(:,T_timepoints(t));
    Theta_color=T_Theta(:,T_timepoints(t));
    Brightness=Zfilt(:,T_timepoints(t));
    for n=1:N
        if Beta_color(n)
            s(n).FaceColor=[(Brightness(n)+1)*0.35 0.5*(Brightness(n)+1) (Brightness(n)+1)*0.35];
            s(n).FaceAlpha=0.7;
        elseif Theta_color(n)
            s(n).FaceColor=[0.5*(Brightness(n)+1) (Brightness(n)+1)*0.35 (Brightness(n)+1)*0.35];
            s(n).FaceAlpha=0.7;
        elseif Alpha_color(n)
            s(n).FaceColor=[(Brightness(n)+1)*0.35 (Brightness(n)+1)*0.35 0.5*(Brightness(n)+1)];
            s(n).FaceAlpha=0.7;
        elseif Delta_color(n)
            s(n).FaceColor=[0.5*(Brightness(n)+1) 0.4*(Brightness(n)+1) 0.25*(Brightness(n)+1)];
            s(n).FaceAlpha=0.7;
        else
            s(n).FaceColor=[0.5+Brightness(n)*.2 0.5+Brightness(n)*.2 0.5+Brightness(n)*.2];
            s(n).FaceAlpha=0.1;
        end
    end
    title(['T = ' num2str(T_timepoints(t)*dt_save,'%3.2f') 's'])
    view(t,20)
    %pause(0.05)
    camlight(h,'headlight')
    frame = getframe(gcf);
%     writeVideo(videoMOMs,frame); %write the image to file
end
close(videoMOMs); %close the file
%SAVE VIDEO



%% FIGURE with sequence of MOMs over top view from Top

T_timepoints=3500:50:6000;

figure('color','w')

for t=1:length(T_timepoints)
    
    subplot_tight(5,ceil(numel(T_timepoints)/5),t,0.001)
    
    hold on
    % First plot a transparent cortex
    cortexpatch=patch(isosurface(scortex,0.1), 'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'none','FaceAlpha', 0.1);
    reducepatch(cortexpatch,0.1);
    isonormals(scortex,cortexpatch);
    
    % Then plot one sphere in each area
    for n=1:N
        s(n)=surf(x+MNI_coord(n,2), y+MNI_coord(n,1),z+MNI_coord(n,3),'FaceColor',[0.9 .9 .9],'EdgeColor','none','FaceAlpha',0.2);
    end
    
    view(90,-90)
    camlight('headlight')
    axis equal
    material dull
    lighting gouraud
    daspect([1 1 1])

    xlim([5 105])
    ylim([5 85])
    zlim([0 80])
    axis off

    Delta_color=T_Delta(:,T_timepoints(t));
    Alpha_color=T_Alpha(:,T_timepoints(t));
    Beta_color=T_Beta(:,T_timepoints(t));
    Theta_color=T_Theta(:,T_timepoints(t));
    Brightness=ones(1,N);
    for n=1:N
        if Beta_color(n)
            s(n).FaceColor=[(Brightness(n)+1)*0.5 0.5*(Brightness(n)+1) (Brightness(n)+1)*0.3];
            s(n).FaceAlpha=0.7;
        elseif Theta_color(n)
            s(n).FaceColor=[0.5*(Brightness(n)+1) (Brightness(n)+1)*0.3 (Brightness(n)+1)*0.3];
            s(n).FaceAlpha=0.7;
        elseif Alpha_color(n)
            s(n).FaceColor=[(Brightness(n)+1)*0.3 (Brightness(n)+1)*0.3 0.5*(Brightness(n)+1)];
            s(n).FaceAlpha=0.7;
        elseif Delta_color(n)
            s(n).FaceColor=[0.5*(Brightness(n)+1) 0.4*(Brightness(n)+1) 0.25*(Brightness(n)+1)];
            s(n).FaceAlpha=0.7;
        else
            s(n).FaceColor=[0.5+Brightness(n)*.2 0.5+Brightness(n)*.2 0.5+Brightness(n)*.2];
            s(n).FaceAlpha=0.1;
        end
    end
    title(['T = ' num2str(T_timepoints(t)*dt_save,'%3.2f') 's'])
    %pause(0.05)
end


