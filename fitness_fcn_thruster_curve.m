function loss = fitness_fcn_thruster_curve(x, Va, Throttle, dt, t_data, T_Data, Q_data)
%FITNESS_FCN_CONTROL_SYSTEM for plotting
%function loss = fitness_fcn_thruster_curve(x, Va, Throttle, dt, t_data, T_Data, Q_data)
%
%Change Optimization Parameters
Thruster_Config.g = 32.5; %[rps/throttle]
%Thruster_Config.k = x(1); %[rate parameter]
Thruster_Config.k1 = x(1); %[rate parameter]
Thruster_Config.k2 = x(2); %[rate parameter]
Thruster_Config.kt = x(3); %[rate parameter]
Thruster_Config.D = 0.1151; %[m] propellor diameter
Thruster_Config.alpha1 = x(4);
Thruster_Config.alpha2 = 0;
Thruster_Config.beta1 = abs(x(5));
Thruster_Config.beta2 = 0;
Thruster_Config.RH_prop = (x(6)>0); %1 for RH, 0 for LH

rho = 1027; %[kg/m^3] Density of seawater

%Run Simulation
[t_n,T_out, Q_out, n_out] = sim_thruster(Va, Throttle , rho, Thruster_Config, dt);

%Interpolate thrust data and calculate mse
T_Data_int = interp1(t_data,T_Data, t_n);
mse_T = mean((T_Data_int-T_out).^2); %Calculate MSE

%Interpolate torque data and calculate mse
Q_Data_int = interp1(t_data,Q_data, t_n);
mse_Q = mean((Q_Data_int-Q_out).^2);

%Calculate loss function
loss = 1*mse_T + 1*mse_Q;
