function mse = fitness_fcn_control_system(x, To, Ta, dt, T_Data);
%FITNESS_FCN_CONTROL_SYSTEM for plotting
%function mse = fitness_fcn_control_system(x,y,Control_Params, Model_Params);
%
%Change Optimization Parameters
k = x(1);

%Run Simulation
[t_n, T_out] = sim_thruster(To, Ta , k, dt);

%Calculate mse between simulation and data 
mse=mean((T_Data-T_out).^2); %Calculate MSE
