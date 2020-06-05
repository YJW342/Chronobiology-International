% Evolutionary Strategies with Evolution
clear all
close all
clc
bdclose('all');
max_iterations=100; % maximum iterations
mu=200;             % numbers of every generation
lambda=100;         % numbers of children
rho=2;             % number of parents to create one children
Mutation=0.05;      % mutation percentage
clear zeta gamma_omg gamma_d
Cost=[];N_zeta_min=-2;N_zeta_max=0;N_omg_min=-25;N_omg_max=-8;N_d_min=-25;N_d_max=0;
% generate random initial guesses in the first generation
N=rand(mu,3);
Zeta=10.^(N(:,1)*(N_zeta_max-N_zeta_min)+N_zeta_min);
Gamma_omg=10.^(N(:,2)*(N_omg_max-N_omg_min)+N_omg_min);
Gamma_d=10.^(N(:,3)*(N_d_max-N_d_min)+N_d_min);
%filename=strcat('A3.mat');
%load(filename);
% sample waveform
Data=[0:1/60:400]';
for i=1:size(Data,1)
    Data(i,2)=1000+500*sin(2*pi*Data(i,1)/24+pi/2)+300*sin(4*pi*Data(i,1)/24+pi)+200*sin(6*pi*Data(i,1)/24+pi/2)+100*sin(8*pi*Data(i,1)/24+3*pi/2)+50*sin(10*pi*Data(i,1)/24)+200*randn;
end
% load('data_fine.mat'); % load .mat data file, the first column is time series and the second column is actigrahy data

T=Data(2,1)-Data(1,1);   % Sampling period
Fs = 1/T;                % Sampling frequency
L = size(Data,1);        % Length of signal
t = Data(:,1)-Data(1,1); % Time vector
X=Data(:,end);Y = fft(X);P2 = abs(Y/L);P = P2(1:L/2+1);P(2:end-1) = 2*P(2:end-1);
f = Fs*(0:(L/2))/L;
x_initial = [0 0 0 0 0 0 0 2*pi/24];
% calculate the cost values of the first generation
for i=1:size(Zeta,1)
    zeta=Zeta(i);gamma_omg=Gamma_omg(i);gamma_d=Gamma_d(i);
    sim('ANF_3rd.mdl');
    Y=fft(y_ANF(:,2));P2 = abs(Y/L);P1 = P2(1:L/2+1);P1(2:end-1) = 2*P1(2:end-1);
    [M1,N1]=min(abs(f-1/24));[M2,N2]=min(abs(f-0.0289));NN=N1-N2;
    J_harm=trapz((P1(1:NN)-P(1:NN)).^2)+trapz((P1(N1-NN:N1+NN)-P(N1-NN:N1+NN)).^2);
    J_noise=trapz(P1(NN+1:N1-NN-1).^2)+trapz(P1(N1+NN+1:end).^2);
    Cost(i)=J_harm+J_noise;
end
for i=1:max_iterations
    Combination=nchoosek([1:50],rho);Labels=randperm(size(Combination,1));
    for j=1:lambda
        Mutation_children=rand(3,1);
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
        Zeta=[Zeta;Feature_Children(1)];
        Gamma_omg=[Gamma_omg;Feature_Children(2)];
        Gamma_d=[Gamma_d;Feature_Children(3)];
        zeta=Zeta(end);gamma_omg=Gamma_omg(end);gamma_d=Gamma_d(end);
        sim('ANF_3rd.mdl');
        % sim('ANF_3rd.mdl');
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
zeta=mean(Zeta)
gamma_omg=mean(Gamma_omg)
gamma_d=mean(Gamma_d)
    


    