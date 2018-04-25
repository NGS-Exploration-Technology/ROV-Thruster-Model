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
n_samples = 3000;
Throttle = ones(1,n_samples);
Va = zeros(size(Throttle));
dt = 0.001;

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
[t, Fx, Mx] = process_raw_thruster_data();

t_Data = t;
T_Data = Fx;
Q_Data = Mx;


%define Fitness Function
%f = @(x)fitness_fcn_thruster_curve(x, To, Ta, dt, t_Data, Data);
f = @(x)fitness_fcn_thruster_curve(x, Va, Throttle, dt, t_Data, T_Data, Q_Data);

[x fval exitflag opt_output] = fminsearch(f,[ 1 1 1 1 1 1], options)

%simulate and plot result
Thruster_Params.g = 32.5; %[rps/throttle]
Thruster_Params.k = x(1); %[rate parameter]
Thruster_Params.D = 0.1151; %[m] propellor diameter
Thruster_Params.alpha1 = x(2);
Thruster_Params.alpha2 = x(3);
Thruster_Params.beta1 = x(4);
Thruster_Params.beta2 = x(5);
rho = 1027;

%Run Simulation
[t_fit,T_fit, Q_fit, n_fit] = sim_thruster(Va, Throttle , rho, Thruster_Params, dt);

%Plot results
figure;
subplot(3,1,1);
plot(t_Data, T_Data, t_fit, T_fit, '--k'); 
ylabel('Force [N]');
legend('Experimental', 'Model');
axis([0 3 0 100]);
grid on;
subplot(3,1,2);
plot(t_Data,Q_Data, t_fit, Q_fit, '--k');
xlabel('Time [s]'); ylabel('Torque [Nm]');
axis([0 3 -1 20]);
grid on;
subplot(3,1,3);
plot(t_fit, n_fit);
xlabel('Time [s]'); ylabel('Prop Speed [rps]');
axis([0 3 -1 35]);
grid on;