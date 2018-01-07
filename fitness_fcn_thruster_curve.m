function mse = fitness_fcn_control_system(x, To, Ta, dt, t_data, Data);
%FITNESS_FCN_CONTROL_SYSTEM for plotting
%function mse = fitness_fcn_control_system(x,y,Control_Params, Model_Params);
%
%Change Optimization Parameters
k = x(1);

%Run Simulation
[t_n, T_out] = sim_thruster(To, Ta , k, dt);

%Interpolate data for comparison
Data_int = interp1(t_data,Data, t_n);

%Calculate mse between simulation and data 
mse=mean((Data_int-T_out).^2); %Calculate MSE
