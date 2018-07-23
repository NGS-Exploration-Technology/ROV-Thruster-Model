%Set fminsearch options
options.Display= 'final';
options.MaxFunEvals = 5000;
options.MaxIter = [];
options.TolFun = 1E-12;
options.TolX = 1E-12;
options.FunValCheck= [];
options.OutputFcn = [];
options.PlotFcns = {@optimplotfval,@optimplotx};

% Set pso options:
% psooptions = optimoptions('particleswarm');
% psooptions.CreationFcn= @pswcreationuniform;
% psooptions.Display= 'iter';
% psooptions.FunctionTolerance= 1.0000e-12;
% psooptions.HybridFcn= [];
% psooptions.InertiaRange= [0.1000 1.1000];
% psooptions.InitialSwarmMatrix= [];
% psooptions.InitialSwarmSpan= 100000;
% psooptions.MaxIterations= Inf;
% psooptions.MaxStallIterations= 1000;
% psooptions.MaxStallTime= Inf;
% psooptions.MaxTime= Inf;
% psooptions.MinNeighborsFraction= 0.2500;
% psooptions.ObjectiveLimit= -Inf;
% psooptions.OutputFcn=[];
% psooptions.PlotFcn= {[@pswplotbestf]};
% psooptions.SelfAdjustmentWeight= 1.4900;
% psooptions.SocialAdjustmentWeight= 1.4900;
% psooptions.SwarmSize = 30;
% psooptions.UseParallel= false;
% psooptions.UseVectorized= false;

%Data to fit (Newton's law of cooling)
%Ta_Data = 20;
%k_Data = 20;
%t_Data = linspace(0,20,10000);
%Data=Ta_Data+(To-Ta_Data)*exp(-k_Data*t_Data) + 0.1*randn(size(t_Data));

%Data to fit (Charge on a capacitor)
%t = linspace(0,5,10000); %[s]
%C = 2E-6; %[F]
%R = 1E5; %[Ohms]
%Vb = 0.5E6; %[V]
%Q = C*Vb*(1-exp(-t/(R*C))) + 1E-8*randn(size(t));
[t, n, Fx, Mx] = process_raw_thruster_data();

t_Data = t;
n_Data = n;
T_Data = Fx;
Q_Data = Mx;

%Set up model params
n_samples = length(t_Data);
Throttle = 5*ones(1,n_samples);
Va = zeros(size(Throttle));
dt = 0.001;

constants = Generate_Thrust_Curves()
prop_const = optimize_voltage_constant()

%define Fitness Function
%f = @(x)fitness_fcn_thruster_curve(x, To, Ta, dt, t_Data, Data);
f = @(x)fitness_fcn_thruster_curve(x, constants, prop_const, Va, Throttle, dt, t_Data, n_Data, T_Data, Q_Data);

[x fval exitflag opt_output] = fminsearch(f,[700 10 .5*constants(2) 15*constants(6)], options)
% [x fval exitflag opt_output] = particleswarm(f,8,zeros(1,8),500*ones(1,8),psooptions)

% x = [10 10 1 43042.5 .1 .1 .00000000000001 .000000000000001];
% x = [73.9277 .1 .1 1950*73.9277 .1 .1 0 0];
% x = [11.4875 .0001 .437 5152.929 1.38857 4.092 .0001 .000000001];
% x = [1.0894e-13 .0043 70.7035 4676.2 2.2955 179.7885 5.5639e-06 8.8988e-08];
% x = [10 0 79500 5500 0 10 2*constants(2) 2*constants(6)];
% x = [10 0 2500 5500 0 10 2*constants(2) 2*constants(6)];
% x = [600 10 .5*constants(2) 15*constants(6)];

%simulate and plot result
rho = 1027; %[kg/m^3] Density of seawater
Thruster_Config.D = 0.1151; %[m] propellor diameter
Thruster_Config.L = 0.07; %[m] duct length
Thruster_Config.kt = (2*rho*Thruster_Config.L*(pi/4)*(Thruster_Config.D^2))^-1; % Thrust coeff from Fossen paper
Thruster_Config.kn1 = abs(prop_const(1)); %[rate parameter]
Thruster_Config.kv = abs(prop_const(2)); %[rate parameter]
Thruster_Config.kq = abs(x(1)); %[rate parameter]
Thruster_Config.ku1 = abs(x(2)); %[rate parameter]
Thruster_Config.ku2 = Thruster_Config.L^-1; %[rate parameter]
Thruster_Config.cT1 = constants(1);
Thruster_Config.cT2 = constants(2);
Thruster_Config.dT1 = constants(3);
Thruster_Config.dT2 = constants(4);
Thruster_Config.cQ1 = constants(5);
Thruster_Config.cQ2 = constants(6);
Thruster_Config.dQ1 = constants(7);
Thruster_Config.dQ2 = constants(8);
Thruster_Config.alpha2 = abs(x(3));
Thruster_Config.beta2 = abs(x(4));
Thruster_Config.RH_prop = (constants(5)<0); %1 for RH, 0 for LH

%Run Simulation
[t_fit, T_fit, Q_fit, n_fit, u_fit] = sim_thruster(Va, Throttle , rho, Thruster_Config, dt);

%Plot results
figure;
plot(t_Data, T_Data, t_fit, T_fit, '--k');
title('Thrust Profile Model Fit to Data');
xlabel('Time, t, [s]'); ylabel('Thrust, T, [N]');
legend('Experimental', 'Model');
axis([0 1 0 60]);
grid on;

figure;
plot(t_Data, -Q_Data, t_fit, -Q_fit, '--k');
title('Torque Profile Model Fit to Data');
xlabel('Time, t, [s]'); ylabel('Torque [Nm]');
legend('Experimental', 'Model');
axis([0 1 0 2]);
grid on;

figure;
plot(t_Data, n_Data, t_fit, n_fit, '--k');
title('Propeller Speed Model Fit to Data')
xlabel('Time, t, [s]'); ylabel('Prop Speed, n, [rpm]');
legend('Experimental', 'Model');
axis([0 1 0 2500]);
grid on;

figure;
plot(t_fit,u_fit,'--k')
xlabel('Time [s]'); ylabel('Flow Speed [m/s]');
grid on;