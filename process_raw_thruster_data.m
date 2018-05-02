function [t, Fx, Mx] = process_raw_thruster_data()
%PROCESS_RAW_THRUSTER_DATA process and filter raw thruster data
%function [t, Fx, Mx] = process_raw_thruster_data()

A = csvread('201803071252.csv');

t_raw = A(:,2);
t = t_raw(1:5000);

Fx_raw = -A(:,14);
Fx_movmean = movmean(Fx_raw,20);
Fx = Fx_movmean(6442:11441);

Mx_raw = -A(:,17); %RH Thruster Data
%Mx_raw = A(:,17); %LH Thruster Data
Mx_movmean = movmean(Mx_raw,20);
Mx = Mx_movmean(6442:11441);

%subplot(2,1,1);
%plot(t, Fx);

%subplot(2,1,2);
%plot(t, Mx);