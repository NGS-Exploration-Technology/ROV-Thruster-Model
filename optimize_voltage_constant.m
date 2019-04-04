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

% Initialize some parameters:
lenp5 = 4635:7250;
dv1 = -.8475;
dv2 = .9254;


% Read in data:
Ap5 = csvread('201807201712_Air_RPM_Data.csv');
tp5 = Ap5(1:length(lenp5),2);
np5 = (477.43*Ap5(lenp5,21)+5.3255)*2*pi/60;
vp5 = [linspace(dv2,5,41).';5*ones(length(tp5)-632,1);linspace(5,dv2,41).';dv2*ones(550,1)];
% t = A(1:1701,2);
% n = -477.43*A(6800:8500,21)+5.3255;
% v = -5*ones(length(t),1);
% t = A(1:1721,2);
% n = 477.43*A(3180:4900,21)+5.3255;
% v = 3*ones(length(t),1);
dt = tp5(2)-tp5(1);

yyaxis left
plot(tp5,np5)
yyaxis right
plot(tp5,vp5)
grid

% fminsearch operation:
f = @(x)fitness_fcn_rpm_voltage(x,tp5,np5,vp5,dt);

[x,fval,exitflag,opt_output] = fminsearch(f,[10 700], options)

% Set constants:
kn = x(1);
kv2 = x(2);

% Calculate model:
n_fit = zeros(length(tp5),1);
for i = 2:length(n_fit)
    d_n = -kn*n_fit(i-1)+kv2*(vp5(i-1)-dv2);
    n_fit(i) = n_fit(i-1)+d_n*dt;
end

figure
plot(tp5,np5,tp5,n_fit,'--k')
legend('Data','Model')
xlabel('Time, t, [sec]')
ylabel('Rotor Speed, n, [rpm]')
grid

end