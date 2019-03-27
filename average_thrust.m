%% Average Thrust Histories
%% Initialize
clc; clear; close all;

%% Run Routine
lpFilt = designfilt('lowpassiir','FilterOrder',4, ...
         'PassbandFrequency',15,'PassbandRipple',0.01, ...
         'StopbandAttenuation', 120, 'SampleRate',1e3);

Set1 = csvread('Fwd_Rvrs_Set1.csv');
Set2 = csvread('Fwd_Rvrs_Set2.csv');
Set3 = csvread('Fwd_Rvrs_Set3.csv');
Set4 = csvread('Fwd_Rvrs_Set4.csv');
Set5 = csvread('Fwd_Rvrs_Set5.csv');

len1 = 770:5270;
len2 = 2440:6940;
len3 = 810:5310;
len4 = 200:4700;
len5 = 910:5410;

t1 = Set1(len1,2);
t2 = Set2(len2,2);
t3 = Set3(len3,2);
t4 = Set4(len4,2);
t5 = Set5(len5,2);

n1 = (477.43*Set1(len1,21)+5.3255)*2*pi/60;
n2 = (477.43*Set2(len2,21)+5.3255)*2*pi/60;
n3 = (477.43*Set3(len3,21)+5.3255)*2*pi/60;
n4 = (477.43*Set4(len4,21)+5.3255)*2*pi/60;
n5 = (477.43*Set5(len5,21)+5.3255)*2*pi/60;

F1 = Set1(len1,16)*-1;
F2 = Set2(len2,16)*-1;
F3 = Set3(len3,16)*-1;
F4 = Set4(len4,16)*-1;
F5 = Set5(len5,16)*-1;

% figure
% yyaxis left
% plot(t1,F1,t2,F2,t3,F3,t4,F4,t5,F5)
% axis([.5 5.5 -20 20])
% yyaxis right
% plot(t5,n5)
% grid

Fmat = [F1.';F2.';F3.';F4.';F5.'];
Favg = mean(Fmat);
nmat = [n1.';n2.';n3.';n4.';n5.'];
navg = mean(nmat);
Flpf = filter(lpFilt,Favg);

figure
plot(1:4501,Favg)
grid