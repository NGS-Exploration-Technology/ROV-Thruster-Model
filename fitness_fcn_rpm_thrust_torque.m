function loss = fitness_fcn_rpm_thrust_torque(x,n_data,T_data,Q_data)
% Initialize optimization parameters:
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
n2 = n_data.*abs(n_data);

% Delete problematic data:
n2Q = n2;
n2Q(7) = [];
Q_data(7) = [];

% Calculate deadband piecewise quadratic function:
alpha1 = [cT1*(n2(n2<=dT1)-dT1);0*n2(n2<dT2 & n2>dT1);cT2*(n2(n2>=dT2)-dT2)];
beta1 = [cQ1*(n2Q(n2Q<=dQ1)-dQ1);0*n2Q(n2Q<dQ2 & n2Q>dQ1);cQ2*(n2Q(n2Q>=dQ2)-dQ2)];

% Calculate Thrust/Torque:
T = rho*D^4*alpha1; % [N] thrust
Q = rho*D^5*beta1; % [Nm] torque

% Calculate cost function:
loss = sqrt(mean((T-T_data).^2)+mean((Q-Q_data).^2));
end