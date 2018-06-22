function [t, Fx, Mx] = process_raw_thruster_data()
%PROCESS_RAW_THRUSTER_DATA process and filter raw thruster data
%function [t, Fx, Mx] = process_raw_thruster_data()

A = csvread('201803071252.csv');

t_raw = A(:,2);
t = t_raw(1:5000);

Fx_raw = -A(:,14);
%Fx_movmean = movmean(Fx_raw,20);
%Fx = Fx_movmean(6442:11441);

Mx_raw = -A(:,17); %RH Thruster Data
%Mx_raw = A(:,17); %LH Thruster Data
%Mx_movmean = movmean(Mx_raw,20);
%Mx = Mx_movmean(6442:11441);

%subplot(2,1,1);
%plot(t, Fx);

%subplot(2,1,2);
%plot(t, Mx);

%bsf Moment
bsFilt = designfilt('bandstopfir', 'FilterOrder',40, ...
         'CutoffFrequency1',70,'CutoffFrequency2',130, ...
         'StopbandAttenuation', 120, 'SampleRate',1000);
     
Mx_bsf = filter(bsFilt,Mx_raw);
     
%lpf Moment
lpFilt = designfilt('lowpassiir','FilterOrder',20, ...
         'PassbandFrequency',25,'PassbandRipple',0.2, ...
         'StopbandAttenuation', 80, 'SampleRate',1e3);

Mx_bsf_lpf = filter(lpFilt,Mx_bsf);

Mx = Mx_bsf_lpf(6496:11495);

%figure;
%subplot(3,1,1);
%plot(t_raw,Mx_raw); ylabel('[Nm]'); grid on;
%subplot(3,1,2);
%plot(t_raw,Mx_bsf); ylabel('[Nm]'); grid on;
%subplot(3,1,3); 
%plot(t_raw,Mx_bsf_lpf); xlabel('time [s]'); ylabel('[Nm]'); grid on;

%bsf Force
bsFilt = designfilt('bandstopfir', 'FilterOrder',40, ...
         'CutoffFrequency1',70,'CutoffFrequency2',130, ...
         'StopbandAttenuation', 120, 'SampleRate',1000);
     
Fx_bsf = filter(bsFilt,Fx_raw);
     
%lpf Force
lpFilt = designfilt('lowpassiir','FilterOrder',20, ...
         'PassbandFrequency',20,'PassbandRipple',0.2, ...
         'StopbandAttenuation', 120, 'SampleRate',1e3);

Fx_bsf_lpf = filter(lpFilt,Fx_bsf);

Fx = Fx_bsf_lpf(6496:11495);


%figure;
%subplot(3,1,1);
%plot(t_raw,Fx_raw); ylabel('[N]'); grid on;
%subplot(3,1,2);
%plot(t_raw,Fx_bsf); ylabel('[N]'); grid on;
%subplot(3,1,3); 
%plot(t_raw,Fx_bsf_lpf); xlabel('time [s]'); ylabel('[N]'); grid on;

%y = fft(Fx_raw);
%fs = 1e3; %Hz
%n = length(Fx_raw);          % number of samples
%f = (0:n-1)*(fs/n);     % frequency range
%power = abs(y).^2/n;    % power of the DFT

%figure; plot(f,power);
%keyboard;
