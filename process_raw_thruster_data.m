function [t, n, Fx, Mx] = process_raw_thruster_data()
%PROCESS_RAW_THRUSTER_DATA process and filter raw thruster data
%function [t, Fx, Mx] = process_raw_thruster_data()

A = csvread('201807161552_Final_With_Tach.csv');

start_sample = 2533;
stop_sample = 8000;
sample_length = stop_sample-start_sample + 1;

t_raw = A(:,2);
t = t_raw(1:sample_length);

n_raw = 477.43*A(:,21)+5.3255; % Linear regression estimate of rpm from Tach Voltage

lpFilt = designfilt('lowpassiir','FilterOrder',4, ...
         'PassbandFrequency',2,'PassbandRipple',0.01, ...
         'StopbandAttenuation', 120, 'SampleRate',1e3);

lpFilt_reduced = designfilt('lowpassiir','FilterOrder',4, ...
         'PassbandFrequency',15,'PassbandRipple',0.01, ...
         'StopbandAttenuation', 120, 'SampleRate',1e3);
     
% n_lpf = filter(lpFilt,n_raw);

% n = n_lpf(2600:8000);

n_median_filt = median_filter(n_raw, 10);

%n = n_raw(2520:7920);
n = n_median_filt(start_sample:stop_sample);

% n = movmean(n_raw(2500:12000),70);

% figure
% subplot(2,1,1);
% plot(t_raw(1:5481),n_raw(2520:8000)); ylabel('[rpm]'); grid on;
% subplot(2,1,2);
% plot(t,n)
% xlabel('time, [s]'),ylabel('[rpm]')
%grid

% keyboard

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

Mx_mf = median_filter(Mx_raw, 10);
Mx_mf_lpf = filter(lpFilt_reduced,Mx_mf);

Mx_lpf = filter(lpFilt,Mx_raw);
%Mx = Mx_lpf(start_sample:stop_sample);
Mx = Mx_mf_lpf(start_sample:stop_sample);
%plot(t,Mx_raw(start_sample:stop_sample),t,Mx,t,Mx_mf_lpf(start_sample:stop_sample));

% figure;
% subplot(2,1,1);
% plot(t_raw(1:5501),Mx_raw(2500:8000)); ylabel('[Nm]'); grid on;
% subplot(2,1,2); 
% plot(t,Mx); xlabel('time [s]'); ylabel('[Nm]'); grid on;

Fx_lpf = filter(lpFilt,Fx_raw);

Fx_mf = median_filter(Fx_raw, 10);
Fx_mf_lpf = filter(lpFilt_reduced,Fx_mf);
%Fx = Fx_lpf(start_sample:stop_sample);

%plot(t,Fx_raw(start_sample:stop_sample),t,Fx,t,Fx_mf_lpf(start_sample:stop_sample));

Fx = Fx_mf_lpf(start_sample:stop_sample);

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