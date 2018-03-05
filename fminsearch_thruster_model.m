%Set fminsearch options
options.Display= 'final';
options.MaxFunEvals = [];
options.MaxIter = [];
options.TolFun = 1E-12;
options.TolX = 1E-12;
options.FunValCheck= [];
options.OutputFcn = [];
options.PlotFcns = {[@optimplotfval],[ @optimplotx]};

%Set up model params
n_samples = 4000;
Throttle = ones(1,n_samples);
Va = zeros(size(Throttle));
dt = 0.001;

%Data to fit (Newton's law of cooling)
%Ta_Data = 20;
%k_Data = 20;
%t_Data = linspace(0,20,10000);
%Data=Ta_Data+(To-Ta_Data)*exp(-k_Data*t_Data) + 0.1*randn(size(t_Data));

%Data to fit (Charge on a capacitor)
t = linspace(0,5,10000); %[s]
C = 2E-6; %[F]
R = 1E5; %[Ohms]
Vb = 0.5E6; %[V]
Q = C*Vb*(1-exp(-t/(R*C))) + 1E-8*randn(size(t));
t_Data = t;
T_Data = Q;
Q_Data = zeros(size(T_Data));


%define Fitness Function
%f = @(x)fitness_fcn_thruster_curve(x, To, Ta, dt, t_Data, Data);
f = @(x)fitness_fcn_thruster_curve(x, Va, Throttle, dt, t_Data, T_Data, Q_Data);

[x fval exitflag opt_output] = fminsearch(f,[ 4.9875 1 0 0 0], options)

%simulate and plot result
Thruster_Params.g = 1; %[rps/throttle]
Thruster_Params.k = x(1); %[rate parameter]
Thruster_Params.D = 0.05; %[m] propellor diameter
Thruster_Params.alpha1 = x(2);
Thruster_Params.alpha2 = x(3);
Thruster_Params.beta1 = x(4);
Thruster_Params.beta2 = x(5);
rho = 1027;

%Run Simulation
[t_fit,T_fit, Q_fit] = sim_thruster(Va, Throttle , rho, Thruster_Params, dt);

%Plot results
figure;
subplot(2,1,1);
plot(t_Data, T_Data, t_fit, T_fit, '.k'); 
title('Thrust versus Time'); xlabel('Time [s]'); ylabel('Force [N]');
subplot(2,1,2);
plot(t_Data,Q_Data, t_fit, Q_fit, '.k');
title('Torque versus Time'); ylabel('Time [s]'); ylabel('Torque [Nm]');