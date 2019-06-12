% Initialize
clear; close all;

% Import data
lenstepL = 560:10560;
lensineL = 550:10550;
lenstep1 = 400:10400;
lensine1 = 500:10500;
lenstep2kf1 = 270:10270;
lensine2kf1 = 430:10430;
t = 0:.04:10;
load('CLLookupStep.mat')
ndLstep = nd;
TdLstep = Td;
uLstep = u;
yLstep = y;
StepL = csvread('CLLookupStep.csv');
tLstep = StepL(1:length(lenstepL),2);
TLstep = -StepL(lenstepL,16);
load('CLLookupSine.mat')
ndLsine = nd;
TdLsine = Td;
uLsine = u;
yLsine = y;
SineL = csvread('CLLookupSine.csv');
tLsine = SineL(1:length(lensineL),2);
TLsine = -SineL(lensineL,16);
load('CL1stateStep.mat')
nd1step = nd;
nhat1step = nhat(1:(end-1));
P1step = P(1:(end-1));
Td1step = Td;
u1step = u;
y1step = y;
Step1 = csvread('CL1stateStep.csv');
t1step = Step1(1:length(lenstep1),2);
T1step = -Step1(lenstep1,16);
load('CL1stateSine.mat')
nd1sine = nd;
nhat1sine = nhat(1:(end-1));
P1sine = P(1:(end-1));
Td1sine = Td;
u1sine = u;
y1sine = y;
Sine1 = csvread('CL1stateSine.csv');
t1sine = Sine1(1:length(lensine1),2);
T1sine = -Sine1(lensine1,16);
load('CL2stateKF1Step.mat')
nd2KF1step = nd;
nhat2KF1step = nhat(1:(end-1));
vhat2KF1step = vhat(1:(end-1));
vahat2KF1step = vahat(1:(end-1));
P2KF1step = P(:,:,1:(end-1));
Td2KF1step = Td;
u2KF1step = u;
y2KF1step = y;
Step2KF1 = csvread('CL2stateKF1Step.csv');
t2KF1step = Step2KF1(1:length(lenstep2kf1),2);
T2KF1step = -Step2KF1(lenstep2kf1,16);
load('CL2stateKF1Sine.mat')
nd2KF1sine = nd;
nhat2KF1sine = nhat(1:(end-1));
vhat2KF1sine = vhat(1:(end-1));
vahat2KF1sine = vahat(1:(end-1));
P2KF1sine = P(:,:,1:(end-1));
Td2KF1sine = Td;
u2KF1sine = u;
y2KF1sine = y;
Sine2KF1 = csvread('CL2stateKF1Sine.csv');
t2KF1sine = Sine2KF1(1:length(lensine2kf1),2);
T2KF1sine = -Sine2KF1(lensine2kf1,16);

% Plot
figure
subplot(2,3,1)
plot(tLstep,TLstep,t,TdLstep,'--k')
ylabel('Step Thrust [N]')
title('Lookup Table Thrust Results')
grid
subplot(2,3,4)
plot(tLsine,TLsine,t,TdLsine,'--k')
xlabel('Time [s]'),ylabel('Sinusoidal Thrust [N]')
grid
subplot(2,3,2)
plot(t1step,T1step,t,Td1step,'--k')
title('Single-State CL Thrust Results')
grid
subplot(2,3,5)
plot(t1sine,T1sine,t,Td1sine,'--k')
xlabel('Time [s]')
grid
subplot(2,3,3)
plot(t2KF1step,T2KF1step,t,Td2KF1step,'--k')
legend('Data','Setpoint')
title('Two-State EKF Model 1 CL Thrust Results')
grid
subplot(2,3,6)
plot(t2KF1sine,T2KF1sine,t,Td2KF1sine,'--k')
xlabel('Time [s]')
grid

figure
subplot(2,3,1)
plot(t,sign(ndLstep).*yLstep,t,ndLstep,'--k')
ylabel('Step Angular Velocity [rad/s]')
title('Lookup Table Angular Velocity Results')
grid
subplot(2,3,4)
plot(t,sign(ndLsine).*yLsine,t,ndLsine,'--k')
xlabel('Time [s]'),ylabel('Sinusoidal Angular Velocity [rad/s]')
grid
subplot(2,3,2)
plot(t,sign(nhat1step).*y1step,t,nhat1step,t,nd1step,'--k',t,nhat1step+3*sqrt(P1step),'--g',t,nhat1step-3*sqrt(P1step),'--g')
title('Single-State CL Angular Velocity Results')
grid
subplot(2,3,5)
plot(t,sign(nhat1sine).*y1sine,t,nhat1sine,t,nd1sine,'--k',t,nhat1sine+3*sqrt(P1sine),'--g',t,nhat1sine-3*sqrt(P1sine),'--g')
xlabel('Time [s]')
grid
subplot(2,3,3)
plot(t,sign(nhat2KF1step).*(y2KF1step(1,:).'),t,nhat2KF1step,t,nd2KF1step,'--k',t,nhat2KF1step+3*sqrt(squeeze(P2KF1step(1,1,:))),'--g',t,nhat2KF1step-3*sqrt(squeeze(P2KF1step(1,1,:))),'--g')
legend('Measurement','Estimate','Setpoint','3\sigma Bounds')
title('Two-State EKF Model 1 CL Angular Velocity Results')
grid
subplot(2,3,6)
plot(t,sign(nhat2KF1sine).*(y2KF1sine(1,:).'),t,nhat2KF1sine,t,nd2KF1sine,'--k',t,nhat2KF1sine+3*sqrt(squeeze(P2KF1sine(1,1,:))),'--g',t,nhat2KF1sine-3*sqrt(squeeze(P2KF1sine(1,1,:))),'--g')
xlabel('Time [s]')
grid

figure
subplot(2,2,1)
plot(t,vhat2KF1step,t,vhat2KF1step+3*sqrt(squeeze(P2KF1step(2,2,:))),'--g',t,vhat2KF1step-3*sqrt(squeeze(P2KF1step(2,2,:))),'--g')
ylabel('Step Axial Velocity [rad/s]')
title('Lookup Table Angular Velocity Results')
grid
subplot(2,1,2)
plot(t,y2KF1step(1,:),t,sign(ndLsine).*yLsine,t,ndLsine,'--k')
xlabel('Time [s]'),ylabel('Sinusoidal Angular Velocity [rad/s]')
grid