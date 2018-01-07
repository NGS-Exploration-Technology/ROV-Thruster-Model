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
k_Data = 1;
t_Data = linspace(0,20,10000);
Data=Ta(1)+(To-Ta(1))*exp(-k_Data*t_Data) + 0.01*randn(size(t_Data));

t_fit = [0:dt:floor(max(t_Data)/dt)*dt];

%define Fitness Function
f = @(x)fitness_fcn_thruster_curve(x, To, Ta, dt, t_Data, Data);

[x fval exitflag opt_output] = fminsearch(f,[1], options)

%simulate and plot result
k = x(1);
[t_fit, T_fit] = sim_thruster(To, Ta , k, dt);


%Plot results
figure;
plot(t_Data, Data, t_fit, T_fit, '.k');
