Thrust_to_Throttle_Table = csvread('Thrust_to_Throttle_Table.csv');
Table_Thrusts = Thrust_to_Throttle_Table(:,1);
Table_Throttles = Thrust_to_Throttle_Table(:,2);

Thrust = 16; %[N]

Throttle = 5*thrust_to_throttle(Thrust, Table_Thrusts, Table_Throttles);

%Set up model params
n_samples = 10000;
Throttle_vector = Throttle*ones(1,n_samples);
Va = zeros(size(Throttle_vector));
dt = 0.001; %[s]
rho = 1027; %[kg/m^3] Density of seawater


Thruster_Params = setup_thrusters_standard_ROV_full_vector_rev6;
Thruster_Config = Thruster_Params(1); %Parameters for one thruster

[t_fit, T_fit, Q_fit, n_fit] = sim_thruster_throttle(Throttle_vector, rho, Thruster_Config, dt);

%Plot results
figure;
subplot(3,1,1);
plot(t_fit, T_fit, '--k'); 
ylabel('Force [N]');
grid on;

subplot(3,1,2);
plot(t_fit, Q_fit, '--k');
ylabel('Torque [Nm]');
grid on;

subplot(3,1,3);
plot(t_fit, n_fit, '--k');
ylabel('Prop Speed [rpm]');
grid on;
