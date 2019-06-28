%% Step 1: Setup

% Initialize
clear; close all;

% Preallocate indexing parameters
lenstepL = 560:10560;
lensineL = 550:10550;
lenrampL = 560:10560;
lenstep1 = 400:10400;
lensine1 = 500:10500;
lenramp1 = 610:10610;
lenstep2kf1 = 270:10270;
lensine2kf1 = 430:10430;
lenramp2kf1 = 450:10450;
lenstep2nlo1 = 380:10380;
lensine2nlo1 = 660:10660;
lenramp2nlo1 = 300:10300;
lenstep2kfnlo1 = 390:10390;
lensine2kfnlo1 = 430:10430;
lenramp2kfnlo1 = 600:10600;
lenstep2kf2 = 550:10550;
lensine2kf2 = 420:10420;
lenramp2kf2 = 340:10340;
lenstep2nlo2 = 600:10600;
lensine2nlo2 = 500:10500;
lenramp2nlo2 = 520:10520;
lenstep2kfnlo2 = 390:10390;
lensine2kfnlo2 = 470:10470;
lenramp2kfnlo2 = 300:10300;
t_cntr = 0:.04:10;
t_data = 0:.001:10;
Tdstep = [-10*ones((length(t_data)-1)/2,1);10*ones((length(t_data)+1)/2,1)];
Tdsine = 5*cos(1.88*t_data.')-5;
Tdramp = -[(linspace(0,10,2000)).';10*ones(1000,1);(linspace(10,-10,4000)).';-10*ones(1001,1);(linspace(-10,0,2000)).'];
steplim1 = 1200;
steplim2 = 6200;
ramplim1 = 600;
ramplim2 = 4800;
ramplim3 = 6000;

%% Step 2a: Import Step Response Data

% Single State data
load('CLLookupStep.mat')
ndLstep = nd;
uLstep = u;
yLstep = y;
StepL = csvread('CLLookupStep.csv');
TLstep = -StepL(lenstepL,16);
load('CL1stateStep.mat')
nd1step = nd;
nhat1step = nhat(1:(end-1));
P1step = P(1:(end-1));
u1step = u;
y1step = y;
Step1 = csvread('CL1stateStep.csv');
T1step = -Step1(lenstep1,16);

% Hydrodynamic Linear data
load('CL2stateKF1Step.mat')
nd2KF1step = nd;
nhat2KF1step = nhat(1:(end-1));
vhat2KF1step = vhat(1:(end-1));
vahat2KF1step = vahat(1:(end-1));
P2KF1step = P(:,:,1:(end-1));
u2KF1step = u;
y2KF1step = y;
Step2KF1 = csvread('CL2stateKF1Step.csv');
T2KF1step = -Step2KF1(lenstep2kf1,16);
load('CL2stateNLO1Step.mat')
nd2NLO1step = nd;
nhat2NLO1step = nhat(1:(end-1));
vhat2NLO1step = vhat(1:(end-1));
vahat2NLO1step = vahat(1:(end-1));
P2NLO1step = P(1:(end-1));
u2NLO1step = u;
y2NLO1step = y;
Step2NLO1 = csvread('CL2stateNLO1Step.csv');
T2NLO1step = -Step2NLO1(lenstep2nlo1,16);
load('CL2stateKFNLO1Step.mat')
nd2KFNLO1step = nd;
nhat2KFNLO1step = nhat(1:(end-1));
vhat2KFNLO1step = vhat(1:(end-1));
vahat2KFNLO1step = vahat(1:(end-1));
P2KFNLO1step = P(:,:,1:(end-1));
u2KFNLO1step = u;
y2KFNLO1step = y;
Step2KFNLO1 = csvread('CL2stateKFNLO1Step.csv');
T2KFNLO1step = -Step2KFNLO1(lenstep2kfnlo1,16);

% Hydrodynamic Quadratic data
load('CL2stateKF2Step.mat')
nd2KF2step = nd;
nhat2KF2step = nhat(1:(end-1));
vhat2KF2step = vhat(1:(end-1));
vahat2KF2step = vahat(1:(end-1));
P2KF2step = P(:,:,1:(end-1));
u2KF2step = u;
y2KF2step = y;
Step2KF2 = csvread('CL2stateKF2Step.csv');
T2KF2step = -Step2KF2(lenstep2kf2,16);
load('CL2stateNLO2Step.mat')
nd2NLO2step = nd;
nhat2NLO2step = nhat(1:(end-1));
vhat2NLO2step = vhat(1:(end-1));
vahat2NLO2step = vahat(1:(end-1));
P2NLO2step = P(1:(end-1));
u2NLO2step = u;
y2NLO2step = y;
Step2NLO2 = csvread('CL2stateNLO2Step.csv');
T2NLO2step = -Step2NLO2(lenstep2nlo2,16);
load('CL2stateKFNLO2Step.mat')
nd2KFNLO2step = nd;
nhat2KFNLO2step = nhat(1:(end-1));
vhat2KFNLO2step = vhat(1:(end-1));
vahat2KFNLO2step = vahat(1:(end-1));
P2KFNLO2step = P(:,:,1:(end-1));
u2KFNLO2step = u;
y2KFNLO2step = y;
Step2KFNLO2 = csvread('CL2stateKFNLO2Step.csv');
T2KFNLO2step = -Step2KFNLO2(lenstep2kfnlo2,16);

% Calculate filtered data
T1f = [NaN(steplim1,1);T1step((steplim1+1):((length(T1step)-1)/2));NaN(steplim2-((length(T1step)-1)/2),1);T1step((steplim2+1):end)];
TLf = [NaN(steplim1,1);TLstep((steplim1+1):((length(TLstep)-1)/2));NaN(steplim2-((length(TLstep)-1)/2),1);TLstep((steplim2+1):end)];
TKF1f = [NaN(steplim1,1);T2KF1step((steplim1+1):((length(T2KF1step)-1)/2));NaN(steplim2-((length(T2KF1step)-1)/2),1);T2KF1step((steplim2+1):end)];
TKF2f = [NaN(steplim1,1);T2KF2step((steplim1+1):((length(T2KF2step)-1)/2));NaN(steplim2-((length(T2KF2step)-1)/2),1);T2KF2step((steplim2+1):end)];
TNLO1f = [NaN(steplim1,1);T2NLO1step((steplim1+1):((length(T2NLO1step)-1)/2));NaN(steplim2-((length(T2NLO1step)-1)/2),1);T2NLO1step((steplim2+1):end)];
TNLO2f = [NaN(steplim1,1);T2NLO2step((steplim1+1):((length(T2NLO2step)-1)/2));NaN(steplim2-((length(T2NLO2step)-1)/2),1);T2NLO2step((steplim2+1):end)];
TKFNLO1f = [NaN(steplim1,1);T2KFNLO1step((steplim1+1):((length(T2KFNLO1step)-1)/2));NaN(steplim2-((length(T2KFNLO1step)-1)/2),1);T2KFNLO1step((steplim2+1):end)];
TKFNLO2f = [NaN(steplim1,1);T2KFNLO2step((steplim1+1):((length(T2KFNLO2step)-1)/2));NaN(steplim2-((length(T2KFNLO2step)-1)/2),1);T2KFNLO2step((steplim2+1):end)];

% Calculate errors from step data
rmseL = sqrt(nanmean((TLf-Tdstep).^2));
rmse1 = sqrt(nanmean((T1f-Tdstep).^2))/rmseL;
rmse2KF1 = sqrt(nanmean((TKF1f-Tdstep).^2))/rmseL;
rmse2NLO1 = sqrt(nanmean((TNLO1f-Tdstep).^2))/rmseL;
rmse2KFNLO1 = sqrt(nanmean((TKFNLO1f-Tdstep).^2))/rmseL;
rmse2KF2 = sqrt(nanmean((TKF2f-Tdstep).^2))/rmseL;
rmse2NLO2 = sqrt(nanmean((TNLO2f-Tdstep).^2))/rmseL;
rmse2KFNLO2 = sqrt(nanmean((TKFNLO2f-Tdstep).^2))/rmseL;
steptbl = [rmse1 rmse2KF1 rmse2KF2;NaN rmse2NLO1 rmse2NLO2;NaN rmse2KFNLO1 rmse2KFNLO2];

%% Step 2b: Plot Step Response Data

% Plot thrust results
figure
subplot(3,3,1)
plot(t_data,T1step,t_data,T1f,t_data,Tdstep,'--k')
ylim([-15 15])
ylabel('EKF Thrust [N]')
title('Single-State Thrust Results')
grid
subplot(3,3,[4 7])
plot(t_data,TLstep,t_data,TLf,t_data,Tdstep,'--k')
ylim([-15 15])
xlabel('Time [s]'),ylabel('Lookup Table Thrust [N]')
grid
subplot(3,3,2)
plot(t_data,T2KF1step,t_data,TKF1f,t_data,Tdstep,'--k')
ylim([-15 15])
title('Hydrodynamic Linear Model Thrust Results')
grid
subplot(3,3,5)
plot(t_data,T2NLO1step,t_data,TNLO1f,t_data,Tdstep,'--k')
ylim([-15 15])
ylabel('NLO Thrust [N]')
grid
subplot(3,3,8)
plot(t_data,T2KFNLO1step,t_data,TKFNLO1f,t_data,Tdstep,'--k')
ylim([-15 15])
xlabel('Time [s]'),ylabel('EKFNLO Thrust [N]')
grid
subplot(3,3,3)
plot(t_data,T2KF2step,t_data,TKF2f,t_data,Tdstep,'--k')
ylim([-15 15])
title('Hydrodynamic Quadratic Model Thrust Results')
grid
subplot(3,3,6)
plot(t_data,T2NLO2step,t_data,TNLO2f,t_data,Tdstep,'--k')
ylim([-15 15])
grid
subplot(3,3,9)
plot(t_data,T2KFNLO2step,t_data,TKFNLO2f,t_data,Tdstep,'--k')
ylim([-15 15])
xlabel('Time [s]')
legend({'Experiment','Corrected','Desired'},'Location','Southeast')
grid
suptitle('Thrust Step Response Results')

% Plot angular velocity results
figure
subplot(3,3,1)
plot(t_cntr,sign(nhat1step).*y1step,t_cntr,nhat1step,t_cntr,nd1step,'--k',t_cntr,nhat1step+3*sqrt(P1step),'--g',t_cntr,nhat1step-3*sqrt(P1step),'--g')
ylim([-120 120])
title('Single-State Angular Velocity Results')
ylabel('EKF Angular Velocity [rad/s]')
grid
subplot(3,3,[4 7])
plot(t_cntr,sign(ndLstep).*yLstep,t_cntr,ndLstep,'--k')
ylim([-120 120])
xlabel('Time [s]'),ylabel('Lookup Table Angular Velocity [rad/s]')
grid
subplot(3,3,2)
plot(t_cntr,sign(nhat2KF1step).*(y2KF1step(1,:).'),t_cntr,nhat2KF1step,t_cntr,nd2KF1step,'--k',t_cntr,nhat2KF1step+3*sqrt(squeeze(P2KF1step(1,1,:))),'--g',t_cntr,nhat2KF1step-3*sqrt(squeeze(P2KF1step(1,1,:))),'--g')
ylim([-120 120])
title('Hydrodynamic Linear Model Angular Velocity Results')
grid
subplot(3,3,5)
plot(t_cntr,sign(nhat2NLO1step).*(y2NLO1step(1,:).'),t_cntr,nhat2NLO1step,t_cntr,nd2NLO1step,'--k',t_cntr,nhat2NLO1step+3*sqrt(P2NLO1step),'--g',t_cntr,nhat2NLO1step-3*sqrt(P2NLO1step),'--g')
ylim([-120 120])
ylabel('NLO Angular Velocity [rad/s]')
grid
subplot(3,3,8)
plot(t_cntr,sign(nhat2KFNLO1step).*(y2KFNLO1step(1,:).'),t_cntr,nhat2KFNLO1step,t_cntr,nd2KFNLO1step,'--k',t_cntr,nhat2KFNLO1step+3*sqrt(squeeze(P2KFNLO1step(1,1,:))),'--g',t_cntr,nhat2KFNLO1step-3*sqrt(squeeze(P2KFNLO1step(1,1,:))),'--g')
ylim([-120 120])
xlabel('Time [s]'),ylabel('EKFNLO Angular Velocity [rad/s]')
grid
subplot(3,3,3)
plot(t_cntr,sign(nhat2KF2step).*(y2KF2step(1,:).'),t_cntr,nhat2KF2step,t_cntr,nd2KF2step,'--k',t_cntr,nhat2KF2step+3*sqrt(squeeze(P2KF2step(1,1,:))),'--g',t_cntr,nhat2KF2step-3*sqrt(squeeze(P2KF2step(1,1,:))),'--g')
ylim([-120 120])
title('Hydrodynamic Quadratic Model Angular Velocity Results')
grid
subplot(3,3,6)
plot(t_cntr,sign(nhat2NLO2step).*(y2NLO2step(1,:).'),t_cntr,nhat2NLO2step,t_cntr,nd2NLO2step,'--k',t_cntr,nhat2NLO2step+3*sqrt(P2NLO2step),'--g',t_cntr,nhat2NLO2step-3*sqrt(P2NLO2step),'--g')
ylim([-120 120])
grid
subplot(3,3,9)
plot(t_cntr,sign(nhat2KFNLO2step).*(y2KFNLO2step(1,:).'),t_cntr,nhat2KFNLO2step,t_cntr,nd2KFNLO2step,'--k',t_cntr,nhat2KFNLO2step+3*sqrt(squeeze(P2KFNLO2step(1,1,:))),'--g',t_cntr,nhat2KFNLO2step-3*sqrt(squeeze(P2KFNLO2step(1,1,:))),'--g')
ylim([-120 120])
xlabel('Time [s]')
legend({'Measurement','Estimate','Setpoint','3\sigma Bounds'},'Location','Southeast')
grid
suptitle('Angular Velocity Step Response Results')

% Plot fluid velocity results
figure
subplot(3,4,1)
hold on
plot(t_cntr,vhat2KF1step,'Color',[0.8500 0.3250 0.0980])
plot(t_cntr,vhat2KF1step+3*sqrt(squeeze(P2KF1step(2,2,:))),'--g',t_cntr,vhat2KF1step-3*sqrt(squeeze(P2KF1step(2,2,:))),'--g')
ylim([-.3 .3])
ylabel('EKF Velocity [m/s]')
title('Hydrodynamic Linear Model Axial Velocity Results')
grid
hold off
subplot(3,4,2)
plot(t_cntr,y2KF1step(2,:),t_cntr,vahat2KF1step,t_cntr,vahat2KF1step+3*sqrt(squeeze(P2KF1step(3,3,:))),'--g',t_cntr,vahat2KF1step-3*sqrt(squeeze(P2KF1step(3,3,:))),'--g')
ylim([-.3 .3])
title('Hydrodynamic Linear Model Ambient Velocity Results')
grid
subplot(3,4,3)
hold on
plot(t_cntr,vhat2KF2step,'Color',[0.8500 0.3250 0.0980])
plot(t_cntr,vhat2KF2step+3*sqrt(squeeze(P2KF2step(2,2,:))),'--g',t_cntr,vhat2KF2step-3*sqrt(squeeze(P2KF2step(2,2,:))),'--g')
ylim([-.3 .3])
title('Hydrodynamic Quadratic Model Axial Velocity Results')
grid
hold off
subplot(3,4,4)
plot(t_cntr,y2KF2step(2,:),t_cntr,vahat2KF2step,t_cntr,vahat2KF2step+3*sqrt(squeeze(P2KF2step(3,3,:))),'--g',t_cntr,vahat2KF2step-3*sqrt(squeeze(P2KF2step(3,3,:))),'--g')
ylim([-.3 .3])
title('Hydrodynamic Quadratic Model Ambient Velocity Results')
grid
subplot(3,4,5)
plot(t_cntr,vhat2NLO1step,'Color',[0.8500 0.3250 0.0980])
ylim([-.3 .3])
ylabel('NLO Velocity [m/s]')
grid
subplot(3,4,6)
plot(t_cntr,y2NLO1step(2,:),t_cntr,vahat2NLO1step)
ylim([-.3 .3])
grid
subplot(3,4,7)
plot(t_cntr,vhat2NLO2step,'Color',[0.8500 0.3250 0.0980])
ylim([-.3 .3])
grid
subplot(3,4,8)
plot(t_cntr,y2NLO2step(2,:),t_cntr,vahat2NLO2step)
ylim([-.3 .3])
grid
subplot(3,4,9)
hold on
plot(t_cntr,vhat2KFNLO1step,'Color',[0.8500 0.3250 0.0980])
plot(t_cntr,vhat2KFNLO1step+3*sqrt(squeeze(P2KFNLO1step(2,2,:))),'--g',t_cntr,vhat2KFNLO1step-3*sqrt(squeeze(P2KFNLO1step(2,2,:))),'--g')
ylim([-.3 .3])
xlabel('Time [s]'),ylabel('EKFNLO Velocity [m/s]')
grid
hold off
subplot(3,4,10)
plot(t_cntr,y2KFNLO1step(2,:),t_cntr,vahat2KFNLO1step,t_cntr,vahat2KFNLO1step+3*sqrt(squeeze(P2KFNLO1step(3,3,:))),'--g',t_cntr,vahat2KFNLO1step-3*sqrt(squeeze(P2KFNLO1step(3,3,:))),'--g')
ylim([-.3 .3])
xlabel('Time [s]')
grid
subplot(3,4,11)
hold on
plot(t_cntr,vhat2KFNLO2step,'Color',[0.8500 0.3250 0.0980])
plot(t_cntr,vhat2KFNLO2step+3*sqrt(squeeze(P2KFNLO2step(2,2,:))),'--g',t_cntr,vhat2KFNLO2step-3*sqrt(squeeze(P2KFNLO2step(2,2,:))),'--g')
ylim([-.3 .3])
xlabel('Time [s]')
grid
hold off
subplot(3,4,12)
plot(t_cntr,y2KFNLO2step(2,:),t_cntr,vahat2KFNLO2step,t_cntr,vahat2KFNLO2step+3*sqrt(squeeze(P2KFNLO2step(3,3,:))),'--g',t_cntr,vahat2KFNLO2step-3*sqrt(squeeze(P2KFNLO2step(3,3,:))),'--g')
ylim([-.3 .3])
xlabel('Time [s]')
legend('Measurement','Estimate','3\sigma Bounds')
grid
suptitle('Fluid Velocity Step Response Results')

% Heatmap of errors
figure
heatmap({'Single-State','Hydrodynamic Linear','Hydrodynamic Quadratic'},{'EKF','NLO','EKFNLO'},steptbl);
xlabel('Dynamic Model'),ylabel('Estimator Type')
title('Thrust Step Response Normalized RMS Tracking Errors')

%% Step 3a: Import Sinusoidal Response Data

% Single State data
load('CLLookupSine.mat')
ndLsine = nd;
uLsine = u;
yLsine = y;
SineL = csvread('CLLookupSine.csv');
TLsine = -SineL(lensineL,16);
load('CL1stateSine.mat')
nd1sine = nd;
nhat1sine = nhat(1:(end-1));
P1sine = P(1:(end-1));
u1sine = u;
y1sine = y;
Sine1 = csvread('CL1stateSine.csv');
T1sine = -Sine1(lensine1,16);

% Hydrodynamic Linear data
load('CL2stateKF1Sine.mat')
nd2KF1sine = nd;
nhat2KF1sine = nhat(1:(end-1));
vhat2KF1sine = vhat(1:(end-1));
vahat2KF1sine = vahat(1:(end-1));
P2KF1sine = P(:,:,1:(end-1));
u2KF1sine = u;
y2KF1sine = y;
Sine2KF1 = csvread('CL2stateKF1Sine.csv');
T2KF1sine = -Sine2KF1(lensine2kf1,16);
load('CL2stateNLO1Sine.mat')
nd2NLO1sine = nd;
nhat2NLO1sine = nhat(1:(end-1));
vhat2NLO1sine = vhat(1:(end-1));
vahat2NLO1sine = vahat(1:(end-1));
P2NLO1sine = P(1:(end-1));
u2NLO1sine = u;
y2NLO1sine = y;
Sine2NLO1 = csvread('CL2stateNLO1Sine.csv');
T2NLO1sine = -Sine2NLO1(lensine2nlo1,16);
load('CL2stateKFNLO1Sine.mat')
nd2KFNLO1sine = nd;
nhat2KFNLO1sine = nhat(1:(end-1));
vhat2KFNLO1sine = vhat(1:(end-1));
vahat2KFNLO1sine = vahat(1:(end-1));
P2KFNLO1sine = P(:,:,1:(end-1));
u2KFNLO1sine = u;
y2KFNLO1sine = y;
Sine2KFNLO1 = csvread('CL2stateKFNLO1Sine.csv');
T2KFNLO1sine = -Sine2KFNLO1(lensine2kfnlo1,16);

% Hydrodynamic Quadratic data
load('CL2stateKF2Sine.mat')
nd2KF2sine = nd;
nhat2KF2sine = nhat(1:(end-1));
vhat2KF2sine = vhat(1:(end-1));
vahat2KF2sine = vahat(1:(end-1));
P2KF2sine = P(:,:,1:(end-1));
u2KF2sine = u;
y2KF2sine = y;
Sine2KF2 = csvread('CL2stateKF2Sine.csv');
T2KF2sine = -Sine2KF2(lensine2kf2,16);
load('CL2stateNLO2Sine.mat')
nd2NLO2sine = nd;
nhat2NLO2sine = nhat(1:(end-1));
vhat2NLO2sine = vhat(1:(end-1));
vahat2NLO2sine = vahat(1:(end-1));
P2NLO2sine = P(1:(end-1));
u2NLO2sine = u;
y2NLO2sine = y;
Sine2NLO2 = csvread('CL2stateNLO2Sine.csv');
T2NLO2sine = -Sine2NLO2(lensine2nlo2,16);
load('CL2stateKFNLO2Sine.mat')
nd2KFNLO2sine = nd;
nhat2KFNLO2sine = nhat(1:(end-1));
vhat2KFNLO2sine = vhat(1:(end-1));
vahat2KFNLO2sine = vahat(1:(end-1));
P2KFNLO2sine = P(:,:,1:(end-1));
u2KFNLO2sine = u;
y2KFNLO2sine = y;
Sine2KFNLO2 = csvread('CL2stateKFNLO2Sine.csv');
T2KFNLO2sine = -Sine2KFNLO2(lensine2kfnlo2,16);

% Calculate errors from sine data
rmseL = sqrt(mean((TLsine-Tdsine).^2));
rmse1 = sqrt(mean((T1sine-Tdsine).^2))/rmseL;
rmse2KF1 = sqrt(mean((T2KF1sine-Tdsine).^2))/rmseL;
rmse2NLO1 = sqrt(mean((T2NLO1sine-Tdsine).^2))/rmseL;
rmse2KFNLO1 = sqrt(mean((T2KFNLO1sine-Tdsine).^2))/rmseL;
rmse2KF2 = sqrt(mean((T2KF2sine-Tdsine).^2))/rmseL;
rmse2NLO2 = sqrt(mean((T2NLO2sine-Tdsine).^2))/rmseL;
rmse2KFNLO2 = sqrt(mean((T2KFNLO2sine-Tdsine).^2))/rmseL;
sinetbl = [rmse1 rmse2KF1 rmse2KF2;NaN rmse2NLO1 rmse2NLO2;NaN rmse2KFNLO1 rmse2KFNLO2];

%% Step 3b: Plot Sinusoidal Response Data

% Plot thrust results
figure
subplot(3,3,1)
plot(t_data,T1sine,t_data,Tdsine,'--k')
ylim([-15 0])
ylabel('EKF Thrust [N]')
title('Single-State Thrust Results')
grid
subplot(3,3,[4 7])
plot(t_data,TLsine,t_data,Tdsine,'--k')
ylim([-15 0])
xlabel('Time [s]'),ylabel('Lookup Table Thrust [N]')
grid
subplot(3,3,2)
plot(t_data,T2KF1sine,t_data,Tdsine,'--k')
ylim([-15 0])
title('Hydrodynamic Linear Model Thrust Results')
grid
subplot(3,3,5)
plot(t_data,T2NLO1sine,t_data,Tdsine,'--k')
ylim([-15 0])
ylabel('NLO Thrust [N]')
grid
subplot(3,3,8)
plot(t_data,T2KFNLO1sine,t_data,Tdsine,'--k')
ylim([-15 0])
xlabel('Time [s]'),ylabel('EKFNLO Thrust [N]')
grid
subplot(3,3,3)
plot(t_data,T2KF2sine,t_data,Tdsine,'--k')
ylim([-15 0])
title('Hydrodynamic Quadratic Model Thrust Results')
grid
subplot(3,3,6)
plot(t_data,T2NLO2sine,t_data,Tdsine,'--k')
ylim([-15 0])
grid
subplot(3,3,9)
plot(t_data,T2KFNLO2sine,t_data,Tdsine,'--k')
ylim([-15 0])
xlabel('Time [s]')
legend({'Experiment','Desired'},'Location','Southeast')
grid
suptitle('Thrust Sinusoidal Response Results')

% Plot angular velocity results
figure
subplot(3,3,1)
plot(t_cntr,sign(nhat1sine).*y1sine,t_cntr,nhat1sine,t_cntr,nd1sine,'--k',t_cntr,nhat1sine+3*sqrt(P1sine),'--g',t_cntr,nhat1sine-3*sqrt(P1sine),'--g')
ylim([-120 0])
title('Single-State Angular Velocity Results')
ylabel('EKF Angular Velocity [rad/s]')
grid
subplot(3,3,[4 7])
plot(t_cntr,sign(ndLsine).*yLsine,t_cntr,ndLsine,'--k')
ylim([-120 0])
xlabel('Time [s]'),ylabel('Lookup Table Angular Velocity [rad/s]')
grid
subplot(3,3,2)
plot(t_cntr,sign(nhat2KF1sine).*(y2KF1sine(1,:).'),t_cntr,nhat2KF1sine,t_cntr,nd2KF1sine,'--k',t_cntr,nhat2KF1sine+3*sqrt(squeeze(P2KF1sine(1,1,:))),'--g',t_cntr,nhat2KF1sine-3*sqrt(squeeze(P2KF1sine(1,1,:))),'--g')
ylim([-120 0])
title('Hydrodynamic Linear Model Angular Velocity Results')
grid
subplot(3,3,5)
plot(t_cntr,sign(nhat2NLO1sine).*(y2NLO1sine(1,:).'),t_cntr,nhat2NLO1sine,t_cntr,nd2NLO1sine,'--k',t_cntr,nhat2NLO1sine+3*sqrt(P2NLO1sine),'--g',t_cntr,nhat2NLO1sine-3*sqrt(P2NLO1sine),'--g')
ylim([-120 0])
ylabel('NLO Angular Velocity [rad/s]')
grid
subplot(3,3,8)
plot(t_cntr,sign(nhat2KFNLO1sine).*(y2KFNLO1sine(1,:).'),t_cntr,nhat2KFNLO1sine,t_cntr,nd2KFNLO1sine,'--k',t_cntr,nhat2KFNLO1sine+3*sqrt(squeeze(P2KFNLO1sine(1,1,:))),'--g',t_cntr,nhat2KFNLO1sine-3*sqrt(squeeze(P2KFNLO1sine(1,1,:))),'--g')
ylim([-120 0])
xlabel('Time [s]'),ylabel('EKFNLO Angular Velocity [rad/s]')
grid
subplot(3,3,3)
plot(t_cntr,sign(nhat2KF2sine).*(y2KF2sine(1,:).'),t_cntr,nhat2KF2sine,t_cntr,nd2KF2sine,'--k',t_cntr,nhat2KF2sine+3*sqrt(squeeze(P2KF2sine(1,1,:))),'--g',t_cntr,nhat2KF2sine-3*sqrt(squeeze(P2KF2sine(1,1,:))),'--g')
ylim([-120 0])
title('Hydrodynamic Quadratic Model Angular Velocity Results')
grid
subplot(3,3,6)
plot(t_cntr,sign(nhat2NLO2sine).*(y2NLO2sine(1,:).'),t_cntr,nhat2NLO2sine,t_cntr,nd2NLO2sine,'--k',t_cntr,nhat2NLO2sine+3*sqrt(P2NLO2sine),'--g',t_cntr,nhat2NLO2sine-3*sqrt(P2NLO2sine),'--g')
ylim([-120 0])
grid
subplot(3,3,9)
plot(t_cntr,sign(nhat2KFNLO2sine).*(y2KFNLO2sine(1,:).'),t_cntr,nhat2KFNLO2sine,t_cntr,nd2KFNLO2sine,'--k',t_cntr,nhat2KFNLO2sine+3*sqrt(squeeze(P2KFNLO2sine(1,1,:))),'--g',t_cntr,nhat2KFNLO2sine-3*sqrt(squeeze(P2KFNLO2sine(1,1,:))),'--g')
ylim([-120 0])
xlabel('Time [s]')
legend('Measurement','Estimate','Setpoint','3\sigma Bounds')
grid
suptitle('Angular Velocity Sinusoidal Response Results')

% Plot fluid velocity results
figure
subplot(3,4,1)
hold on
plot(t_cntr,vhat2KF1sine,'Color',[0.8500 0.3250 0.0980])
plot(t_cntr,vhat2KF1sine+3*sqrt(squeeze(P2KF1sine(2,2,:))),'--g',t_cntr,vhat2KF1sine-3*sqrt(squeeze(P2KF1sine(2,2,:))),'--g')
ylim([-.3 .3])
ylabel('EKF Velocity [m/s]')
title('Hydrodynamic Linear Model Axial Velocity Results')
grid
hold off
subplot(3,4,2)
plot(t_cntr,y2KF1sine(2,:),t_cntr,vahat2KF1sine,t_cntr,vahat2KF1sine+3*sqrt(squeeze(P2KF1sine(3,3,:))),'--g',t_cntr,vahat2KF1sine-3*sqrt(squeeze(P2KF1sine(3,3,:))),'--g')
ylim([-.3 .3])
title('Hydrodynamic Linear Model Ambient Velocity Results')
grid
subplot(3,4,3)
hold on
plot(t_cntr,vhat2KF2sine,'Color',[0.8500 0.3250 0.0980])
plot(t_cntr,vhat2KF2sine+3*sqrt(squeeze(P2KF2sine(2,2,:))),'--g',t_cntr,vhat2KF2sine-3*sqrt(squeeze(P2KF2sine(2,2,:))),'--g')
ylim([-.3 .3])
title('Hydrodynamic Quadratic Model Axial Velocity Results')
grid
hold off
subplot(3,4,4)
plot(t_cntr,y2KF2sine(2,:),t_cntr,vahat2KF2sine,t_cntr,vahat2KF2sine+3*sqrt(squeeze(P2KF2sine(3,3,:))),'--g',t_cntr,vahat2KF2sine-3*sqrt(squeeze(P2KF2sine(3,3,:))),'--g')
ylim([-.3 .3])
title('Hydrodynamic Quadratic Model Ambient Velocity Results')
grid
subplot(3,4,5)
plot(t_cntr,vhat2NLO1sine,'Color',[0.8500 0.3250 0.0980])
ylim([-.3 .3])
ylabel('NLO Velocity [m/s]')
grid
subplot(3,4,6)
plot(t_cntr,y2NLO1sine(2,:),t_cntr,vahat2NLO1sine)
ylim([-.3 .3])
grid
subplot(3,4,7)
plot(t_cntr,vhat2NLO2sine,'Color',[0.8500 0.3250 0.0980])
ylim([-.3 .3])
grid
subplot(3,4,8)
plot(t_cntr,y2NLO2sine(2,:),t_cntr,vahat2NLO2sine)
ylim([-.3 .3])
grid
subplot(3,4,9)
hold on
plot(t_cntr,vhat2KFNLO1sine,'Color',[0.8500 0.3250 0.0980])
plot(t_cntr,vhat2KFNLO1sine+3*sqrt(squeeze(P2KFNLO1sine(2,2,:))),'--g',t_cntr,vhat2KFNLO1sine-3*sqrt(squeeze(P2KFNLO1sine(2,2,:))),'--g')
ylim([-.3 .3])
xlabel('Time [s]'),ylabel('EKFNLO Velocity [m/s]')
grid
hold off
subplot(3,4,10)
plot(t_cntr,y2KFNLO1sine(2,:),t_cntr,vahat2KFNLO1sine,t_cntr,vahat2KFNLO1sine+3*sqrt(squeeze(P2KFNLO1sine(3,3,:))),'--g',t_cntr,vahat2KFNLO1sine-3*sqrt(squeeze(P2KFNLO1sine(3,3,:))),'--g')
ylim([-.3 .3])
xlabel('Time [s]')
grid
subplot(3,4,11)
hold on
plot(t_cntr,vhat2KFNLO2sine,'Color',[0.8500 0.3250 0.0980])
plot(t_cntr,vhat2KFNLO2sine+3*sqrt(squeeze(P2KFNLO2sine(2,2,:))),'--g',t_cntr,vhat2KFNLO2sine-3*sqrt(squeeze(P2KFNLO2sine(2,2,:))),'--g')
ylim([-.3 .3])
xlabel('Time [s]')
grid
hold off
subplot(3,4,12)
plot(t_cntr,y2KFNLO2sine(2,:),t_cntr,vahat2KFNLO2sine,t_cntr,vahat2KFNLO2sine+3*sqrt(squeeze(P2KFNLO2sine(3,3,:))),'--g',t_cntr,vahat2KFNLO2sine-3*sqrt(squeeze(P2KFNLO2sine(3,3,:))),'--g')
ylim([-.3 .3])
xlabel('Time [s]')
legend('Measurement','Estimate','3\sigma Bounds')
grid
suptitle('Fluid Velocity Sinusoidal Response Results')

% Heatmap of errors
figure
heatmap({'Single-State','Hydrodynamic Linear','Hydrodynamic Quadratic'},{'EKF','NLO','EKFNLO'},sinetbl);
xlabel('Dynamic Model'),ylabel('Estimator Type')
title('Thrust Sinusoidal Response Normalized RMS Tracking Errors')

%% Step 4a: Import Ramp Response Data

% Single State data
load('CLLookupRamp.mat')
ndLramp = nd;
uLramp = u;
yLramp = y;
RampL = csvread('CLLookupRamp.csv');
TLramp = -RampL(lenrampL,16);
load('CL1stateRamp.mat')
nd1ramp = nd;
nhat1ramp = nhat(1:(end-1));
P1ramp = P(1:(end-1));
u1ramp = u;
y1ramp = y;
Ramp1 = csvread('CL1stateRamp.csv');
T1ramp = -Ramp1(lenramp1,16);

% Hydrodynamic Linear data
load('CL2stateKF1Ramp.mat')
nd2KF1ramp = nd;
nhat2KF1ramp = nhat(1:(end-1));
vhat2KF1ramp = vhat(1:(end-1));
vahat2KF1ramp = vahat(1:(end-1));
P2KF1ramp = P(:,:,1:(end-1));
u2KF1ramp = u;
y2KF1ramp = y;
Ramp2KF1 = csvread('CL2stateKF1Ramp.csv');
T2KF1ramp = -Ramp2KF1(lenramp2kf1,16);
load('CL2stateNLO1Ramp.mat')
nd2NLO1ramp = nd;
nhat2NLO1ramp = nhat(1:(end-1));
vhat2NLO1ramp = vhat(1:(end-1));
vahat2NLO1ramp = vahat(1:(end-1));
P2NLO1ramp = P(1:(end-1));
u2NLO1ramp = u;
y2NLO1ramp = y;
Ramp2NLO1 = csvread('CL2stateNLO1Ramp.csv');
T2NLO1ramp = -Ramp2NLO1(lenramp2nlo1,16);
load('CL2stateKFNLO1Ramp.mat')
nd2KFNLO1ramp = nd;
nhat2KFNLO1ramp = nhat(1:(end-1));
vhat2KFNLO1ramp = vhat(1:(end-1));
vahat2KFNLO1ramp = vahat(1:(end-1));
P2KFNLO1ramp = P(:,:,1:(end-1));
u2KFNLO1ramp = u;
y2KFNLO1ramp = y;
Ramp2KFNLO1 = csvread('CL2stateKFNLO1Ramp.csv');
T2KFNLO1ramp = -Ramp2KFNLO1(lenramp2kfnlo1,16);

% Hydrodynamic Quadratic data
load('CL2stateKF2Ramp.mat')
nd2KF2ramp = nd;
nhat2KF2ramp = nhat(1:(end-1));
vhat2KF2ramp = vhat(1:(end-1));
vahat2KF2ramp = vahat(1:(end-1));
P2KF2ramp = P(:,:,1:(end-1));
u2KF2ramp = u;
y2KF2ramp = y;
Ramp2KF2 = csvread('CL2stateKF2Ramp.csv');
T2KF2ramp = -Ramp2KF2(lenramp2kf2,16);
load('CL2stateNLO2Ramp.mat')
nd2NLO2ramp = nd;
nhat2NLO2ramp = nhat(1:(end-1));
vhat2NLO2ramp = vhat(1:(end-1));
vahat2NLO2ramp = vahat(1:(end-1));
P2NLO2ramp = P(1:(end-1));
u2NLO2ramp = u;
y2NLO2ramp = y;
Ramp2NLO2 = csvread('CL2stateNLO2Ramp.csv');
T2NLO2ramp = -Ramp2NLO2(lenramp2nlo2,16);
load('CL2stateKFNLO2Ramp.mat')
nd2KFNLO2ramp = nd;
nhat2KFNLO2ramp = nhat(1:(end-1));
vhat2KFNLO2ramp = vhat(1:(end-1));
vahat2KFNLO2ramp = vahat(1:(end-1));
P2KFNLO2ramp = P(:,:,1:(end-1));
u2KFNLO2ramp = u;
y2KFNLO2ramp = y;
Ramp2KFNLO2 = csvread('CL2stateKFNLO2Ramp.csv');
T2KFNLO2ramp = -Ramp2KFNLO2(lenramp2kfnlo2,16);

% Calculate filtered data
T1f = [NaN(ramplim1,1);T1ramp((ramplim1+1):ramplim2);NaN(size(T1ramp((ramplim2+1):ramplim3)));T1ramp((ramplim3+1):end)];
TLf = [NaN(ramplim1,1);TLramp((ramplim1+1):ramplim2);NaN(size(TLramp((ramplim2+1):ramplim3)));TLramp((ramplim3+1):end)];
TKF1f = [NaN(ramplim1,1);T2KF1ramp((ramplim1+1):ramplim2);NaN(size(T2KF1ramp((ramplim2+1):ramplim3)));T2KF1ramp((ramplim3+1):end)];
TKF2f = [NaN(ramplim1,1);T2KF2ramp((ramplim1+1):ramplim2);NaN(size(T2KF2ramp((ramplim2+1):ramplim3)));T2KF2ramp((ramplim3+1):end)];
TNLO1f = [NaN(ramplim1,1);T2NLO1ramp((ramplim1+1):ramplim2);NaN(size(T2NLO1ramp((ramplim2+1):ramplim3)));T2NLO1ramp((ramplim3+1):end)];
TNLO2f = [NaN(ramplim1,1);T2NLO2ramp((ramplim1+1):ramplim2);NaN(size(T2NLO2ramp((ramplim2+1):ramplim3)));T2NLO2ramp((ramplim3+1):end)];
TKFNLO1f = [NaN(ramplim1,1);T2KFNLO1ramp((ramplim1+1):ramplim2);NaN(size(T2KFNLO1ramp((ramplim2+1):ramplim3)));T2KFNLO1ramp((ramplim3+1):end)];
TKFNLO2f = [NaN(ramplim1,1);T2KFNLO2ramp((ramplim1+1):ramplim2);NaN(size(T2KFNLO2ramp((ramplim2+1):ramplim3)));T2KFNLO2ramp((ramplim3+1):end)];

% Calculate errors from step data
rmseL = sqrt(nanmean((TLf-Tdramp).^2));
rmse1 = sqrt(nanmean((T1f-Tdramp).^2))/rmseL;
rmse2KF1 = sqrt(nanmean((TKF1f-Tdramp).^2))/rmseL;
rmse2NLO1 = sqrt(nanmean((TNLO1f-Tdramp).^2))/rmseL;
rmse2KFNLO1 = sqrt(nanmean((TKFNLO1f-Tdramp).^2))/rmseL;
rmse2KF2 = sqrt(nanmean((TKF2f-Tdramp).^2))/rmseL;
rmse2NLO2 = sqrt(nanmean((TNLO2f-Tdramp).^2))/rmseL;
rmse2KFNLO2 = sqrt(nanmean((TKFNLO2f-Tdramp).^2))/rmseL;
ramptbl = [rmse1 rmse2KF1 rmse2KF2;NaN rmse2NLO1 rmse2NLO2;NaN rmse2KFNLO1 rmse2KFNLO2];

%% Step 4b: Plot Ramp Response Data

% Plot thrust results
figure
subplot(3,3,1)
plot(t_data,T1ramp,t_data,T1f,t_data,Tdramp,'--k')
ylim([-15 15])
ylabel('EKF Thrust [N]')
title('Single-State Thrust Results')
grid
subplot(3,3,[4 7])
plot(t_data,TLramp,t_data,TLf,t_data,Tdramp,'--k')
ylim([-15 15])
xlabel('Time [s]'),ylabel('Lookup Table Thrust [N]')
grid
subplot(3,3,2)
plot(t_data,T2KF1ramp,t_data,TKF1f,t_data,Tdramp,'--k')
ylim([-15 15])
title('Hydrodynamic Linear Model Thrust Results')
grid
subplot(3,3,5)
plot(t_data,T2NLO1ramp,t_data,TNLO1f,t_data,Tdramp,'--k')
ylim([-15 15])
ylabel('NLO Thrust [N]')
grid
subplot(3,3,8)
plot(t_data,T2KFNLO1ramp,t_data,TKFNLO1f,t_data,Tdramp,'--k')
ylim([-15 15])
xlabel('Time [s]'),ylabel('EKFNLO Thrust [N]')
grid
subplot(3,3,3)
plot(t_data,T2KF2ramp,t_data,TKF2f,t_data,Tdramp,'--k')
ylim([-15 15])
title('Hydrodynamic Quadratic Model Thrust Results')
grid
subplot(3,3,6)
plot(t_data,T2NLO2ramp,t_data,TNLO2f,t_data,Tdramp,'--k')
ylim([-15 15])
grid
subplot(3,3,9)
plot(t_data,T2KFNLO2ramp,t_data,TKFNLO2f,t_data,Tdramp,'--k')
ylim([-15 15])
xlabel('Time [s]')
legend({'Experiment','Corrected','Desired'},'Location','Southeast')
grid
suptitle('Thrust Ramp Response Results')

% Plot angular velocity results
figure
subplot(3,3,1)
plot(t_cntr,sign(nhat1ramp).*y1ramp,t_cntr,nhat1ramp,t_cntr,nd1ramp,'--k',t_cntr,nhat1ramp+3*sqrt(P1ramp),'--g',t_cntr,nhat1ramp-3*sqrt(P1ramp),'--g')
ylim([-120 120])
title('Single-State Angular Velocity Results')
ylabel('EKF Angular Velocity [rad/s]')
grid
subplot(3,3,[4 7])
plot(t_cntr,sign(ndLramp).*yLramp,t_cntr,ndLramp,'--k')
ylim([-120 120])
xlabel('Time [s]'),ylabel('Lookup Table Angular Velocity [rad/s]')
grid
subplot(3,3,2)
plot(t_cntr,sign(nhat2KF1ramp).*(y2KF1ramp(1,:).'),t_cntr,nhat2KF1ramp,t_cntr,nd2KF1ramp,'--k',t_cntr,nhat2KF1ramp+3*sqrt(squeeze(P2KF1ramp(1,1,:))),'--g',t_cntr,nhat2KF1ramp-3*sqrt(squeeze(P2KF1ramp(1,1,:))),'--g')
ylim([-120 120])
title('Hydrodynamic Linear Model Angular Velocity Results')
grid
subplot(3,3,5)
plot(t_cntr,sign(nhat2NLO1ramp).*(y2NLO1ramp(1,:).'),t_cntr,nhat2NLO1ramp,t_cntr,nd2NLO1ramp,'--k',t_cntr,nhat2NLO1ramp+3*sqrt(P2NLO1ramp),'--g',t_cntr,nhat2NLO1ramp-3*sqrt(P2NLO1ramp),'--g')
ylim([-120 120])
ylabel('NLO Angular Velocity [rad/s]')
grid
subplot(3,3,8)
plot(t_cntr,sign(nhat2KFNLO1ramp).*(y2KFNLO1ramp(1,:).'),t_cntr,nhat2KFNLO1ramp,t_cntr,nd2KFNLO1ramp,'--k',t_cntr,nhat2KFNLO1ramp+3*sqrt(squeeze(P2KFNLO1ramp(1,1,:))),'--g',t_cntr,nhat2KFNLO1ramp-3*sqrt(squeeze(P2KFNLO1ramp(1,1,:))),'--g')
ylim([-120 120])
xlabel('Time [s]'),ylabel('EKFNLO Angular Velocity [rad/s]')
grid
subplot(3,3,3)
plot(t_cntr,sign(nhat2KF2ramp).*(y2KF2ramp(1,:).'),t_cntr,nhat2KF2ramp,t_cntr,nd2KF2ramp,'--k',t_cntr,nhat2KF2ramp+3*sqrt(squeeze(P2KF2ramp(1,1,:))),'--g',t_cntr,nhat2KF2ramp-3*sqrt(squeeze(P2KF2ramp(1,1,:))),'--g')
ylim([-120 120])
title('Hydrodynamic Quadratic Model Angular Velocity Results')
grid
subplot(3,3,6)
plot(t_cntr,sign(nhat2NLO2ramp).*(y2NLO2ramp(1,:).'),t_cntr,nhat2NLO2ramp,t_cntr,nd2NLO2ramp,'--k',t_cntr,nhat2NLO2ramp+3*sqrt(P2NLO2ramp),'--g',t_cntr,nhat2NLO2ramp-3*sqrt(P2NLO2ramp),'--g')
ylim([-120 120])
grid
subplot(3,3,9)
plot(t_cntr,sign(nhat2KFNLO2ramp).*(y2KFNLO2ramp(1,:).'),t_cntr,nhat2KFNLO2ramp,t_cntr,nd2KFNLO2ramp,'--k',t_cntr,nhat2KFNLO2ramp+3*sqrt(squeeze(P2KFNLO2ramp(1,1,:))),'--g',t_cntr,nhat2KFNLO2ramp-3*sqrt(squeeze(P2KFNLO2ramp(1,1,:))),'--g')
ylim([-120 120])
xlabel('Time [s]')
legend({'Measurement','Estimate','Setpoint','3\sigma Bounds'},'Location','Southeast')
grid
suptitle('Angular Velocity Ramp Response Results')

% Plot fluid velocity results
figure
subplot(3,4,1)
hold on
plot(t_cntr,vhat2KF1ramp,'Color',[0.8500 0.3250 0.0980])
plot(t_cntr,vhat2KF1ramp+3*sqrt(squeeze(P2KF1ramp(2,2,:))),'--g',t_cntr,vhat2KF1ramp-3*sqrt(squeeze(P2KF1ramp(2,2,:))),'--g')
ylim([-.3 .3])
ylabel('EKF Velocity [m/s]')
title('Hydrodynamic Linear Model Axial Velocity Results')
grid
hold off
subplot(3,4,2)
plot(t_cntr,y2KF1ramp(2,:),t_cntr,vahat2KF1ramp,t_cntr,vahat2KF1ramp+3*sqrt(squeeze(P2KF1ramp(3,3,:))),'--g',t_cntr,vahat2KF1ramp-3*sqrt(squeeze(P2KF1ramp(3,3,:))),'--g')
ylim([-.3 .3])
title('Hydrodynamic Linear Model Ambient Velocity Results')
grid
subplot(3,4,3)
hold on
plot(t_cntr,vhat2KF2ramp,'Color',[0.8500 0.3250 0.0980])
plot(t_cntr,vhat2KF2ramp+3*sqrt(squeeze(P2KF2ramp(2,2,:))),'--g',t_cntr,vhat2KF2ramp-3*sqrt(squeeze(P2KF2ramp(2,2,:))),'--g')
ylim([-.3 .3])
title('Hydrodynamic Quadratic Model Axial Velocity Results')
grid
hold off
subplot(3,4,4)
plot(t_cntr,y2KF2ramp(2,:),t_cntr,vahat2KF2ramp,t_cntr,vahat2KF2ramp+3*sqrt(squeeze(P2KF2ramp(3,3,:))),'--g',t_cntr,vahat2KF2ramp-3*sqrt(squeeze(P2KF2ramp(3,3,:))),'--g')
ylim([-.3 .3])
title('Hydrodynamic Quadratic Model Ambient Velocity Results')
grid
subplot(3,4,5)
plot(t_cntr,vhat2NLO1ramp,'Color',[0.8500 0.3250 0.0980])
ylim([-.3 .3])
ylabel('NLO Velocity [m/s]')
grid
subplot(3,4,6)
plot(t_cntr,y2NLO1ramp(2,:),t_cntr,vahat2NLO1ramp)
ylim([-.3 .3])
grid
subplot(3,4,7)
plot(t_cntr,vhat2NLO2ramp,'Color',[0.8500 0.3250 0.0980])
ylim([-.3 .3])
grid
subplot(3,4,8)
plot(t_cntr,y2NLO2ramp(2,:),t_cntr,vahat2NLO2ramp)
ylim([-.3 .3])
grid
subplot(3,4,9)
hold on
plot(t_cntr,vhat2KFNLO1ramp,'Color',[0.8500 0.3250 0.0980])
plot(t_cntr,vhat2KFNLO1ramp+3*sqrt(squeeze(P2KFNLO1ramp(2,2,:))),'--g',t_cntr,vhat2KFNLO1ramp-3*sqrt(squeeze(P2KFNLO1ramp(2,2,:))),'--g')
ylim([-.3 .3])
xlabel('Time [s]'),ylabel('EKFNLO Velocity [m/s]')
grid
hold off
subplot(3,4,10)
plot(t_cntr,y2KFNLO1ramp(2,:),t_cntr,vahat2KFNLO1ramp,t_cntr,vahat2KFNLO1ramp+3*sqrt(squeeze(P2KFNLO1ramp(3,3,:))),'--g',t_cntr,vahat2KFNLO1ramp-3*sqrt(squeeze(P2KFNLO1ramp(3,3,:))),'--g')
ylim([-.3 .3])
xlabel('Time [s]')
grid
subplot(3,4,11)
hold on
plot(t_cntr,vhat2KFNLO2ramp,'Color',[0.8500 0.3250 0.0980])
plot(t_cntr,vhat2KFNLO2ramp+3*sqrt(squeeze(P2KFNLO2ramp(2,2,:))),'--g',t_cntr,vhat2KFNLO2ramp-3*sqrt(squeeze(P2KFNLO2ramp(2,2,:))),'--g')
ylim([-.3 .3])
xlabel('Time [s]')
grid
hold off
subplot(3,4,12)
plot(t_cntr,y2KFNLO2ramp(2,:),t_cntr,vahat2KFNLO2ramp,t_cntr,vahat2KFNLO2ramp+3*sqrt(squeeze(P2KFNLO2ramp(3,3,:))),'--g',t_cntr,vahat2KFNLO2ramp-3*sqrt(squeeze(P2KFNLO2ramp(3,3,:))),'--g')
ylim([-.3 .3])
xlabel('Time [s]')
legend('Measurement','Estimate','3\sigma Bounds')
grid
suptitle('Fluid Velocity Ramp Response Results')

% Heatmap of errors
figure
heatmap({'Single-State','Hydrodynamic Linear','Hydrodynamic Quadratic'},{'EKF','NLO','EKFNLO'},ramptbl);
xlabel('Dynamic Model'),ylabel('Estimator Type')
title('Thrust Ramp Response Normalized RMS Tracking Errors')