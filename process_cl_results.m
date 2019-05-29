% Initialize
clear;% close all;

% Import data
thruster_config;
cTn = Thruster_Config.cTn;
lensine = 400:10400;
lenstep = 350:10350;
t = 0:.03:10;
Sine = csvread('CLTest1stateSine.csv');
tsine = Sine(1:length(lensine),2);
Tsine = -Sine(lensine,16);
nsine = 49.9723*Sine(lensine,21);
Tdsine = 10*sin(1.88*t.');
ndsine = sqrt(abs(Tdsine/cTn));
Step = csvread('CLTest1stateStep.csv');
tstep = Step(1:length(lenstep),2);
Tstep = -Step(lenstep,16);
nstep = 49.9723*Step(lenstep,21);
Tdstep = [10*ones(length(t)/2,1);-10*ones(length(t)/2,1)];
ndstep = sqrt(abs(Tdstep/cTn));

% Plot
figure
subplot(2,1,1)
plot(tsine,Tsine,t,Tdsine,'--k')
legend('Data','Setpoint')
ylabel('Thrust [N]')
grid
subplot(2,1,2)
plot(tstep,Tstep,t,Tdstep,'--k')
xlabel('Time [s]')
ylabel('Thrust [N]')
grid

figure
subplot(2,1,1)
plot(tsine,nsine,t,ndsine,'--k')
legend('Data','Setpoint')
ylabel('Propeller Velocity [rad/s]')
grid
subplot(2,1,2)
plot(tstep,nstep,t,ndstep,'--k')
xlabel('Time [s]')
ylabel('Propeller Velocity [rad/s]')
grid