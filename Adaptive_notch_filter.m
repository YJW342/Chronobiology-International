%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adaptive notch filter algorithm in our Chronobiology International paper
% 1st-order ANF
clear all
close all
clc
bdclose('all');
% define a sample waveform
Data=[0:1/60:400]';
for i=1:size(Data,1)
    Data(i,2)=1000+500*sin(2*pi*Data(i,1)/24+pi/2)+300*sin(4*pi*Data(i,1)/24+pi)+200*sin(6*pi*Data(i,1)/24+pi/2)+100*sin(8*pi*Data(i,1)/24+3*pi/2)+50*sin(10*pi*Data(i,1)/24)+200*randn;
end
% load('data.mat'); % load .mat data file, the first column is time series and the second column is actigrahy data
t = Data(:,1)-Data(1,1); % Time vector
T = Data(2,1)-Data(1,1); % Time step
Fs = 1/T;                % Sampling frequency
L = size(Data,1);        % Length of signal
X=Data(:,end);Y = fft(X);P2 = abs(Y/L);P = P2(1:L/2+1);P(2:end-1) = 2*P(2:end-1);
f = Fs*(0:(L/2))/L;
x_initial = [0 0 0 2*pi/24]; % initial guess of ANF state
zeta=0.1251; % ANF damping factor
gamma_omg=3.8019*10^(-9); % adaptation rate of frequency
gamma_d=0.0339; % adaptation rate of constant bias
sim('ANF_1st.mdl');
subplot(3,1,1)
plot(Data(:,1),Data(:,2))
hold on
plot(y_ANF(:,1),y_ANF(:,2),'linewidth',4)
subplot(3,1,2)
plot(x(:,1),x(:,2),'linewidth',4)
subplot(3,1,3)
plot(Data(:,1),mod(2*pi*Data(:,1)/24+pi/2,2*pi),'linewidth',4)
hold on
plot(theta(:,1),theta(:,2),'linewidth',4)

% 2nd-order ANF
clc
x_initial = [0 0 0 0 0 2*pi/24]; % initial guess of ANF state
zeta=0.0904; % ANF damping factor
gamma_omg=2.0534*10^(-9); % adaptation rate of frequency
gamma_d=0.0334; % adaptation rate of constant bias
sim('ANF_2nd.mdl');
subplot(3,1,1)
hold on
plot(y_ANF(:,1),y_ANF(:,2),'linewidth',4)
subplot(3,1,2)
hold on
plot(x(:,1),x(:,2),'linewidth',4)
subplot(3,1,3)
hold on
plot(theta(:,1),theta(:,2),'linewidth',4)

% 3rd-order ANF
clc
x_initial = [0 0 0 0 0 0 0 2*pi/24]; % initial guess of ANF state
zeta=0.0750; % ANF damping factor
gamma_omg=4.9957*10^(-10); % adaptation rate of frequency
gamma_d=0.0333; % adaptation rate of constant bias
sim('ANF_3rd.mdl');
subplot(3,1,1)
hold on
plot(y_ANF(:,1),y_ANF(:,2),'linewidth',4)
subplot(3,1,2)
hold on
plot(x(:,1),x(:,2),'linewidth',4)
subplot(3,1,3)
hold on
plot(theta(:,1),theta(:,2),'linewidth',4)

% 4th-order ANF
clc
x_initial = [0 0 0 0 0 0 0 0 0 2*pi/24]; % initial guess of ANF state
zeta=0.0667; % ANF damping factor
gamma_omg=3.5233*10^(-10); % adaptation rate of frequency
gamma_d=0.0329; % adaptation rate of constant bias
sim('ANF_4th.mdl');
subplot(3,1,1)
hold on
plot(y_ANF(:,1),y_ANF(:,2),'linewidth',4)
subplot(3,1,2)
hold on
plot(x(:,1),x(:,2),'linewidth',4)
subplot(3,1,3)
hold on
plot(theta(:,1),theta(:,2),'linewidth',4)

% 5th-order ANF
clc
x_initial = [0 0 0 0 0 0 0 0 0 0 0 2*pi/24]; % initial guess of ANF state
zeta=0.0599; % ANF damping factor
gamma_omg=2.0893*10^(-10); % adaptation rate of frequency
gamma_d=0.0324; % adaptation rate of constant bias
sim('ANF_5th.mdl');
subplot(3,1,1)
hold on
plot(y_ANF(:,1),y_ANF(:,2),'linewidth',4)
grid on
legend('raw data','1st-order ANF','2nd-order ANF','3rd-order ANF','4th-order ANF','5th-order ANF')
ylabel('y(t)')
subplot(3,1,2)
hold on
plot(x(:,1),x(:,2),'linewidth',4)
grid on
legend('1st-order ANF','2nd-order ANF','3rd-order ANF','4th-order ANF','5th-order ANF')
ylabel('x_1(t)')
subplot(3,1,3)
hold on
plot(theta(:,1),theta(:,2),'linewidth',4)
grid on
xlabel('t (hours)')
ylabel('\theta(t)')