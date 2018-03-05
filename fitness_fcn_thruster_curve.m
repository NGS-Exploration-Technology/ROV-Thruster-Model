function loss = fitness_fcn_thruster_curve(x, Va, Throttle, dt, t_data, T_Data, Q_data)
%FITNESS_FCN_CONTROL_SYSTEM for plotting
%function loss = fitness_fcn_thruster_curve(x, Va, Throttle, dt, t_data, T_Data, Q_data)
%
%Change Optimization Parameters
Thruster_Params.g = x(1); %[rps/throttle]
Thruster_Params.k = x(2); %[rate parameter]
Thruster_Params.D = 0.05; %[m] propellor diameter
Thruster_Params.alpha1 = x(3);
Thruster_Params.alpha2 = x(4);
Thruster_Params.beta1 = x(5);
Thruster_Params.beta2 = x(6);

rho = 1027; %[kg/m^3] Density of seawater

%Run Simulation
[t_n,T_out, Q_out] = sim_thruster(Va, Throttle , rho, Thruster_Params, dt);

%Interpolate thrust data and calculate mse
T_Data_int = interp1(t_data,T_Data, t_n);
mse_T = mean((T_Data_int-T_out).^2); %Calculate MSE

%Interpolate torque data and calculate mse
Q_Data_int = interp1(t_data,Q_data, t_n);
mse_Q = mean((Q_Data_int-Q_out).^2);

%Calculate loss function
loss = 1*mse_T + 1*mse_Q;
