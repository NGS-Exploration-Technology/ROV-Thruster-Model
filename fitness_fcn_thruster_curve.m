function loss = fitness_fcn_thruster_curve(x, constants, prop_const, Va, Throttle, dt, t_Data, n_Data, T_Data, Q_Data)
%FITNESS_FCN_CONTROL_SYSTEM for plotting
%function loss = fitness_fcn_thruster_curve(x, Va, Throttle, dt, t_data, T_Data, Q_data)
%
%Change Optimization Parameters
%Thruster_Config.g = 32.5; %[rps/throttle]
rho = 1027; %[kg/m^3] Density of seawater
Thruster_Config.D = 0.1151; %[m] propellor diameter
Thruster_Config.L = 0.075; %[m] duct length
Thruster_Config.kt = (2*rho*.07*(pi/4)*(Thruster_Config.D^2))^-1; % Thrust coeff from Fossen paper
Thruster_Config.kn1 = abs(prop_const(1)); %[rate parameter]
Thruster_Config.kv = abs(prop_const(2)); %[rate parameter]
Thruster_Config.kq = abs(x(1)); %[rate parameter]
Thruster_Config.ku1 = abs(x(2)); %[rate parameter]
Thruster_Config.ku2 = Thruster_Config.L^-1; %[rate parameter]
Thruster_Config.cT1 = constants(1);
Thruster_Config.cT2 = constants(2);
Thruster_Config.dT1 = constants(3);
Thruster_Config.dT2 = constants(4);
Thruster_Config.cQ1 = constants(5);
Thruster_Config.cQ2 = constants(6);
Thruster_Config.dQ1 = constants(7);
Thruster_Config.dQ2 = constants(8);
Thruster_Config.alpha2 = abs(x(3));
Thruster_Config.beta2 = abs(x(4));
Thruster_Config.RH_prop = (constants(5)<0); %1 for RH, 0 for LH

%Run Simulation
[~, T_out, Q_out, n_out, ~] = sim_thruster(Va, Throttle , rho, Thruster_Config, dt);

%Interpolate thrust data and calculate mse
% T_Data_int = interp1(t_Data,T_Data, t_n);
mse_T = mean((T_Data(500:end)-T_out(500:end)').^2); %Calculate MSE
var_T = var(T_Data(500:end)-T_out(500:end)');

%Interpolate torque data and calculate mse
% Q_Data_int = interp1(t_Data,Q_Data, t_n);
mse_Q = mean((Q_Data(500:end)-Q_out(500:end)').^2);
var_Q = var(Q_Data(500:5000)-Q_out(500:5000)');

%Interpolate rpm data and calculate mse
% n_Data_int = interp1(t_Data,n_Data, t_n);
mse_N = mean((n_Data(500:end)-n_out(500:end)').^2);
var_N = var(n_Data(500:end)-n_out(500:end)');

%Calculate loss function
loss = ((40*sqrt(mse_T+var_T)+1000*sqrt(mse_Q+var_Q)+sqrt(mse_N+var_N))/3)+std(sqrt([1600*mse_T 1000000*mse_Q mse_N]));
