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
n = 100;
To = 0;
Ta = ones(1,n);
dt = 0.01;

%Data to fit (Newton's law of cooling)
%Ta_Data = 20;
%k_Data = 20;
%t_Data = linspace(0,20,10000);
%Data=Ta_Data+(To-Ta_Data)*exp(-k_Data*t_Data) + 0.1*randn(size(t_Data));

%Data to fit (Charge on a capacitor)
t = linspace(0,2,10000); %[s]
C = 1E-6; %[F]
R = 1E5; %[Ohms]
Vb = 1E6; %[V]
Q = C*Vb*(1-exp(-t/(R*C))) + 1E-8*randn(size(t));
t_Data = t;
Data = Q;

%Data to fit (Second order system)
%t = linspace(0,2,10000); %[s]
%c1 = 1;
%r1 = 1;
%c2 = 1;
%r2 = -2;
%t_Data = t;
%Data = c1*exp(r1*t)+c2*exp(r2*t);


%define Fitness Function
f = @(x)fitness_fcn_thruster_curve(x, To, Ta, dt, t_Data, Data);

[x fval exitflag opt_output] = fminsearch(f,[2 0], options)

%simulate and plot result
k = x(1);
g = x(2);
I = 0;
Cd = 0;
v = 0;
[t_fit, T_fit] = sim_thruster(To, Ta , k, g, I, Cd, v, dt);

%Plot results
figure;
plot(t_Data, Data, t_fit, T_fit, '.k');
