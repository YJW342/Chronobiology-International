% Evolutionary Strategy with Evolution to optimize the ANF parameter (damping factor 'zeta', adaptation rate for notch frequency 'gamma_omg' and adaptation rate for constant bias 'gamma_d')
% The detail of the Evolutionary Strategy algorithm can be found in "https://www.sft.asso.fr/Local/sft/dir/user-3775/documents/actes/journeessft/metti_5_2012/Lectures&Tutorials-Texts/Text-T2-Ruffio.pdf"
clear all
close all
clc
bdclose('all');
max_iterations=100; % maximum iterations
mu=200;             % numbers of every generation
lambda=100;         % numbers of children
rho=2;              % number of parents to create one children
Mutation=0.05;      % mutation percentage
clear zeta gamma_omg gamma_d
Cost=[]; % define the objective function set
% define the maximum and minimum of the damping factor and adaptation rates
N_zeta_min=-2;N_zeta_max=0;N_omg_min=-25;N_omg_max=-8;N_d_min=-25;N_d_max=0; 
% generate random initial guesses in the first generation
N=rand(mu,3);
Zeta=10.^(N(:,1)*(N_zeta_max-N_zeta_min)+N_zeta_min);
Gamma_omg=10.^(N(:,2)*(N_omg_max-N_omg_min)+N_omg_min);
Gamma_d=10.^(N(:,3)*(N_d_max-N_d_min)+N_d_min);
% load raw data file if you have one, name it as Data in Matlab, the first column is the time sequence and the second column is data value
filename=strcat('A3.mat');
load(filename);
Data=Data(:,[1,end]);
% sample waveform
% Data=[0:1/60:400]';
% for i=1:size(Data,1)
%     Data(i,2)=1000+500*sin(2*pi*Data(i,1)/24+pi/2)+300*sin(4*pi*Data(i,1)/24+pi)+200*sin(6*pi*Data(i,1)/24+pi/2)+100*sin(8*pi*Data(i,1)/24+3*pi/2)+50*sin(10*pi*Data(i,1)/24)+200*randn;
% end
T=Data(2,1)-Data(1,1);   % Sampling period
Fs = 1/T;                % Sampling frequency
L = size(Data,1);        % Length of signal
t = Data(:,1)-Data(1,1); % Time vector
X=Data(:,end);Y = fft(X);P2 = abs(Y/L);P = P2(1:round(L/2)+1);P(2:end-1) = 2*P(2:end-1); % run fft
f = Fs*(0:(L/2))/L;
x_initial = [0 0 0 2*pi/24]; % define the initial guesses of ANF states
% x_initial = [0 0 0 0 0 2*pi/24]; % for ANF_2nd.mdl
% x_initial = [0 0 0 0 0 0 0 2*pi/24]; % for ANF_3rd.mdl
% x_initial = [0 0 0 0 0 0 0 0 0 2*pi/24]; % for ANF_4th.mdl
% x_initial = [0 0 0 0 0 0 0 0 0 0 0 2*pi/24]; % for ANF_5th.mdl
% calculate the cost values of the first generation
for i=1:size(Zeta,1)
    zeta=Zeta(i);gamma_omg=Gamma_omg(i);gamma_d=Gamma_d(i);
    sim('ANF_1st.mdl'); % or ANF_2nd.mdl,ANF_3rd.mdl,ANF_4th.mdl,ANF_5th.mdl
    Y=fft(y_ANF(:,2));P2 = abs(Y/L);P1 = P2(1:L/2+1);P1(2:end-1) = 2*P1(2:end-1);
    [M1,N1]=min(abs(f-1/24));[M2,N2]=min(abs(f-0.0289));NN=N1-N2;
    J_harm=trapz((P1(1:NN)-P(1:NN)).^2)+trapz((P1(N1-NN:N1+NN)-P(N1-NN:N1+NN)).^2);
    J_noise=trapz(P1(NN+1:N1-NN-1).^2)+trapz(P1(N1+NN+1:end).^2);
    Cost(i)=J_harm+J_noise;
end
for i=1:max_iterations
    Combination=nchoosek([1:50],rho);Labels=randperm(size(Combination,1)); % random combine parent generation
    for j=1:lambda
        Mutation_children=rand(3,1); % define mutation rates for three parameter
        Feature_Children=[mean(Zeta(Combination(Labels(j),:))),mean(Gamma_omg(Combination(Labels(j),:))),mean(Gamma_d(Combination(Labels(j),:)))];
        % check the mitation
        if Mutation_children(1)<Mutation;
            Feature_Children(1)=10^(rand*(N_zeta_max-N_zeta_min)+N_zeta_min);
        end
        if Mutation_children(2)<Mutation;
            Feature_Children(2)=10^(rand*(N_omg_max-N_omg_min)+N_omg_min);
        end
        if Mutation_children(3)<Mutation;
            Feature_Children(3)=10^(rand*(N_d_max-N_d_min)+N_d_min);
        end
        % collect all parameters in parent and children generation
        Zeta=[Zeta;Feature_Children(1)];
        Gamma_omg=[Gamma_omg;Feature_Children(2)];
        Gamma_d=[Gamma_d;Feature_Children(3)];
        zeta=Zeta(end);gamma_omg=Gamma_omg(end);gamma_d=Gamma_d(end);
        sim('ANF_1st.mdl'); % or ANF_2nd.mdl,ANF_3rd.mdl,ANF_4th.mdl,ANF_5th.mdl
        Y=fft(y_ANF(:,2));P2 = abs(Y/L);P1 = P2(1:L/2+1);P1(2:end-1) = 2*P1(2:end-1);
        J_harm=trapz((P1(1:NN)-P(1:NN)).^2)+trapz((P1(N1-NN:N1+NN)-P(N1-NN:N1+NN)).^2);
        J_noise=trapz(P1(NN+1:N1-NN-1).^2)+trapz(P1(N1+NN+1:end).^2);
        Cost=[Cost,J_harm+J_noise];
    end
    % eliminate lambda parameters with largest cost values 
    for j=1:lambda
        [M,N]=max(Cost);
        Cost(N)=[];Zeta(N)=[];Gamma_omg(N)=[];Gamma_d(N)=[];
    end
end
% calculate the optimal damping rate and adaptation rates by average value
% in the last generation
zeta=mean(Zeta)
gamma_omg=mean(Gamma_omg)
gamma_d=mean(Gamma_d)
%plot(Zeta)    
 
