function [Zsave, dt_save] = hopf_delays_simu(f,K,MD,SynDelay,sig,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to run simulations of spontaneous whole-brain network activity
% Where each brain area is represented by a dynamical unit with a
% supercritical Hopf bifurcation.
%
% INPUTS
% - a is the bifurcation parameter 
% - f (Hz) is the fundamental frequency of the natural limit cycle (e.g f=40Hz)
% - dt is the temporal resolution: the smaller the dt, the higher the precision 
%   at the expense of a higher cost function.
%
% ** FREE PARAMETERS **
%
% - K is the Global Coupling Weight 
% - MD (seconds) is the mean delay that scales the distance matrix 
% 
%
% All units in Seconds, Meters, Hz
%
% Runs simulations for a preview period t_prev
% Then starts saving Zsave for a time Tmax at dt_save resolution
%
% Example: [Zsave] = Hopf_Delays_Run(f,K,v,sig,C,D,tmax,t_prev,dt_save)
% where [C, D, tmax,t_prev, dt_save] are avariable arguments
%
% Simulation code written by Gustavo Deco 2016 gustavo.deco@upf.edu
% Adapted to complex coordinates by Adrian Ponce
% Adapted to include time-delays by Joana Cabral, January 2017, joana.cabral@psych.ox.ac.uk
% Adapted by Francesca Castaldo, March 2021, francesca.castaldo.20@ucl.ac.uk
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save_data=1;  % Set to 1 to save the simulated data in the end
% When called in a function, it is better to set this to 0.

%% Variable Input Arguments (varargin)
switch numel(varargin)
    case 0 % This mode is useful if only one run is performed
        % If no additional input arguments, then set them here.
        % Define the Structural Network
        
        load SC_90aal_32HCP.mat mat %load your functional connectivity matrices 
        red_mat=mat; 
%       red_mat(red_mat<10)=0; % Removes connections with less than 10
%       detected fibers to speed up
        N=size(red_mat,1);
        %mat(mat<10)=0; 
        C=red_mat/mean(red_mat(ones(N)-eye(N)>0));
        % Such that the mean of all non-diagonal elements is 1.
        
        load SC_90aal_32HCP.mat mat_D %load your distance matrix
        % Distance between areas
        D=mat_D;
        D=D/1000; % Distance matrix in meters
        
        % Define simulation Parameters
        tmax=50; % in seconds
        t_prev=5; % in seconds
        dt_save=2e-3; % Resolution of saved brain activity in seconds
    case 5 % This mode is useful to accelerate simulation speed for loops over the main free parameters
        [C, D, tmax,t_prev, dt_save]=varargin{1:5};
        
end

%% MODEL PARAMETERS 

a=-5;
N=size(C,1); % Number of areas
dt=1e-4;% Resolution of model integration in seconds
iomega=1i*2*pi*f*ones(N,1); % complex frequency in radians/second


kC=K*C*dt; % Scale matrix C with 'K' and 'dt' to avoid doing it at each step
% Note that it was mean(C(:))=1, so now mean(C(:))=K*dt
dsig = sqrt(dt)*sig; % normalize std of noise with dt

% Set a matrix of Delays containing the number of time-steps between areas
% Delays are integer numbers, so depend on the dt considered.
if MD==0
    Delays=zeros(N); % Set all delays to 0.
else
    Delays=(round((D/mean(D(C>0))*MD+SynDelay)/dt));
end
Delays(C==0)=0;

Max_History=max(Delays(:))+1;
% Convert the Delays matrix such that it contains the index of the History that we need to
% retrieve
Delays=Max_History-Delays;

%% Initialization

% Z = History of Z at dt resolution for the length of the longest delay
% History of Z is initialized as complex gaussian noise with mean 0 and std=1.

Z=dt*randn(N,Max_History)+1i*dt*randn(N,Max_History);

% Matrix Z is a NxMax_History matrix that will be continuously updated

% Simulated activity will be saved only every 2 ms
Zsave=zeros(N,tmax/dt_save);
sumz=zeros(N,1);

%% Run simulations
disp(['Now running for K=' num2str(K) ', mean Delay = ' num2str(MD*1e3) 'ms, a=' num2str(a)])

tic % to count the time of each run
nt=0;

for t=dt:dt:t_prev+tmax %t=1:tmax/dt when adding stimulation  % We only start saving after t_prev
    
    Znow=Z(:,end); % The last collumn of Z is 'now'
    
    % Intrinsic dynamics occurring in this dt
    dz= Znow.*(a + iomega - abs(Znow.^2))*dt;
    
    %Coupling Term. This is a way to speed up the code. 
%     for n=1:N
%         sumz(n) = ...
%             sum(kC(n,1:N) .* (Z((1:N) + (N * (Delays(n,1:N) -1))) - Znow(n)));
%     end


    for n=1:N
        sumzn=0; % Intitalize total coupling received into node n
        for p=1:N
            if kC(n,p) % Calculate only input from connected nodes (faster)
                sumzn = sumzn + kC(n,p) * Z(p,Delays(n,p));% - kC(n,p) * Znow(n) ; %It is commented 
% because is better to run this part outside the loop 
            end
        end
        sumz(n)=sumzn - sum(kC(n,:)*Znow(n));
    end
    
    if MD
        % Slide History only if the delays are >0
        Z(:,1:end-1)=Z(:,2:end);
    end
    
    % Update the end of History with evolution in the last time step
    Z(:,end)= Znow + dz + sumz + dsig*randn(N,1) + 1i*dsig*randn(N,1);
    
    % Save dt_save resolution after t_prev
    if ~mod(t,dt_save) && t>t_prev
        nt=nt+1;
        Zsave(:,nt)=Z(:,end);
    end
end

disp(['Finished, lasted ' num2str(round(toc)) ' secs for real ' num2str(t_prev+tmax) ' seconds'])
disp(['Simu speed ' num2str(round(toc/(t_prev+tmax))) ' Realseconds/SimuSecond'])

if save_data
    K_label=num2str(log10(K));
    ind_p=find(K_label=='.');
    if numel(ind_p)
        K_label(ind_p)='p';
    end
    save(['Hopf_Simu_K1E' K_label '_MD' num2str(MD*1e3)],'Zsave','dt_save','K','MD','f','sig','a')
end
