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
n = 200;
To = 0;
Ta = ones(1,n);
dt = 0.01;

%Data to fit
Ta_Data = 20;
k_Data = 20;
t_Data = linspace(0,20,10000);
Data=Ta_Data+(To-Ta_Data)*exp(-k_Data*t_Data) + 0.1*randn(size(t_Data));

t_fit = [0:dt:floor(max(t_Data)/dt)*dt];

%define Fitness Function
f = @(x)fitness_fcn_thruster_curve(x, To, Ta, dt, t_Data, Data);

[x fval exitflag opt_output] = fminsearch(f,[0 0], options)

%simulate and plot result
k = x(1);
g = x(2);
[t_fit, T_fit] = sim_thruster(To, Ta , k, g, dt);


%Plot results
figure;
plot(t_Data, Data, t_fit, T_fit, '.k');
