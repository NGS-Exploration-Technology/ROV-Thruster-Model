function [t, n, Fx, Mx] = process_raw_thruster_data()
%PROCESS_RAW_THRUSTER_DATA process and filter raw thruster data
%function [t, Fx, Mx] = process_raw_thruster_data()

A = csvread('201807161552_Final_With_Tach.csv');

t_raw = A(:,2);
t = t_raw(1:5401);

n_raw = 477.43*A(:,21)+5.3255; % Linear regression estimate of rpm from Tach Voltage

lpFilt = designfilt('lowpassiir','FilterOrder',4, ...
         'PassbandFrequency',2,'PassbandRipple',0.01, ...
         'StopbandAttenuation', 120, 'SampleRate',1e3);

n_lpf = filter(lpFilt,n_raw);

n = n_lpf(2600:8000);

% n = movmean(n_raw(2500:12000),70);

% figure
% subplot(2,1,1);
% plot(t_raw(1:5481),n_raw(2520:8000)); ylabel('[rpm]'); grid on;
% subplot(2,1,2);
% plot(t,n)
% xlabel('time, [s]'),ylabel('[rpm]')
% grid

Fx_raw = A(:,15)*-1; % -1 Multiplier because test was in reverse direction
%Fx_movmean = movmean(Fx_raw,20);
%Fx = Fx_movmean(6442:11441);

Mx_raw = A(:,18)*-1; % -1 Multiplier because test was in reverse direction
%Mx_raw = -A(:,18); %RH Thruster Data
%Mx_raw = A(:,18); %LH Thruster Data
%Mx_movmean = movmean(Mx_raw,20);
%Mx = Mx_movmean(6442:11441);

%subplot(2,1,1);
%plot(t, Fx);

%subplot(2,1,2);
%plot(t, Mx);

Mx_lpf = filter(lpFilt,Mx_raw);
Mx = Mx_lpf(2600:8000);

% figure;
% subplot(2,1,1);
% plot(t_raw(1:5501),Mx_raw(2500:8000)); ylabel('[Nm]'); grid on;
% subplot(2,1,2); 
% plot(t,Mx); xlabel('time [s]'); ylabel('[Nm]'); grid on;

Fx_lpf = filter(lpFilt,Fx_raw);
Fx = Fx_lpf(2600:8000);

% figure;
% subplot(2,1,1);
% plot(t_raw(1:5501),Fx_raw(2500:8000)); ylabel('[N]'); grid on;
% subplot(2,1,2); 
% plot(t,Fx); xlabel('time [s]'); ylabel('[N]'); grid on;

% y = fft(n_raw);
% fs = 1e3; %Hz
% n = length(n_raw);     % number of samples
% f = (0:n-1)*(fs/n);     % frequency range
% power = abs(y).^2/n;    % power of the DFT
% 
% figure; plot(f,power); grid on;
% 
% keyboard;