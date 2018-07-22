function x = optimize_voltage_constant()
% Read in rpm data and fit to first order model

%Set fminsearch options
options.Display= 'final';
options.MaxFunEvals = [];
options.MaxIter = [];
options.TolFun = 1E-12;
options.TolX = 1E-12;
options.FunValCheck= [];
options.OutputFcn = [];
options.PlotFcns = {@optimplotfval,@optimplotx};

% Read in data:
A = csvread('201807201712_Air_RPM_Data.csv');
t = A(4635:6500,2);
n = 477.43*A(4635:6500,21)+5.3255;
v = 5*ones(length(t),1);
dt = t(2)-t(1);

% fminsearch operation:
f = @(x)fitness_fcn_rpm_voltage(x,t,n,v,dt);

[x,fval,exitflag,opt_output] = fminsearch(f,[1 1], options)
x = abs(x);

% Set constants:
k1 = x(1);
kv = x(2);

% Calculate model:
n_fit = zeros(length(t),1);
for i = 2:length(n_fit)
    d_n = -k1*n_fit(i-1)+kv*v(i-1);
    n_fit(i) = n_fit(i-1)+d_n*dt;
end

figure
plot(t,n,t,n_fit,'--k')
legend('Data','Model')
xlabel('Time, t, [sec]')
ylabel('Rotor Speed, n, [rpm]')
grid

end