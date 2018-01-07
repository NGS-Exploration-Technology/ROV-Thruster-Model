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
n = 10;
To = 0;
Ta = ones(1,n);
dt = 0.5;

%Data to fit
k_fit = 1;
t_fit = [0:dt:n*dt];
T_Data=Ta(1)+(To-Ta(1))*exp(-k_fit*t_fit);

%define Fitness Function
f = @(x)fitness_fcn_thruster_curve(x, To, Ta, dt, T_Data);

[x fval exitflag opt_output] = fminsearch(f,[0], options)

%simulate and plot result
k_fit = x(1);
[t_fit, T_fit] = sim_thruster(To, Ta , k_fit, dt);

%calculate fiting function
t = linspace(min(t_fit), max(t_fit), 10000);
k = 1;
Ta_fit = Ta(1);
T=Ta(1)+(To-Ta(1))*exp(-k*t);

%Plot results
figure;
plot(t, T, t_fit, T_fit, '.k');
