function [x, Throttle, n, T_fit, Q_fit] = Generate_Thrust_Curves()
%function [x, Throttle, n, T_fit, Q_fit] = Generate_Thrust_Curves()
%GENERATE_THRUST_CURVES fit Fossen steady state thruster model
%
%Inputs:
%   None (Hard coded to loadcell test data
%
%Outputs
%   x = Steady state model params
%   Throttle = [V] Throttle Voltage
%   n = Thruster rpm
%   T_fit = thrust fit to model
%   Q_fit = reaction torque fit to model
%
%National Geographic Society
%February 4, 2019
%Jordan Boehm (mostly but cleaned up by Eric Berkenpas)

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
A = csvread('Tach_RPM_Thrust_Data.csv');
Throttle = A(:,1);
n = A(:,7);
T = A(:,5);
Q = A(:,6);


% Delete problematic data:
%n2Q = n2;
Throttle(7) = [];
%n2Q(7) = [];
n(7) = [];
T(7) = [];
Q(7) = [];

% fminsearch operation:
f = @(x)fitness_fcn_rpm_thrust_torque(x,n,T,Q);

[x,fval,exitflag,opt_output] = fminsearch(f,[1 1 -1 1 1 1 -1 1], options);

cT1 = x(1);
cT2 = x(2);
dT1 = x(3);
dT2 = x(4);
cQ1 = x(5);
cQ2 = x(6);
dQ1 = x(7);
dQ2 = x(8);

% Initialize constants:
D = 0.1151; %[m] propellor diameter
rho = 1027; %[kg/m^3] Density of seawater

% Precalculate quadratic term:
%n2 = n.*abs(n); -eb

% Delete problematic data:
%n2Q = n2;
%Throttle(7) = [];
%n2Q(7) = [];
%Q(7) = [];
%n(7) = [];

% Precalculate quadratic term:
n2 = n.*abs(n);
n2Q = n2;

% Calculate deadband piecewise quadratic function:
alpha1 = [cT1*(n2(n2<=dT1)-dT1);0*n2(n2<dT2 & n2>dT1);cT2*(n2(n2>=dT2)-dT2)];
beta1 = [cQ1*(n2Q(n2Q<=dQ1)-dQ1);0*n2Q(n2Q<dQ2 & n2Q>dQ1);cQ2*(n2Q(n2Q>=dQ2)-dQ2)];

% Calculate Thrust/Torque:
T_fit = rho*D^4*alpha1; % [N] thrust
Q_fit = rho*D^5*beta1; % [Nm] torque

figure
plot(n,T,'o',n,T_fit,'--')
xlabel('Propeller Speed, n [rpm]')
ylabel('Thrust, T [N]')
title('Fitted Thrust vs. RPM')
legend('Data','Model')
grid;

%n(7) = [];
figure
plot(n,Q,'o',n,Q_fit,'--')
xlabel('Propeller Speed, n [rpm]')
ylabel('Torque, Q [Nm]')
title('Fitted Torque vs. RPM')
legend('Data','Model')
grid

end