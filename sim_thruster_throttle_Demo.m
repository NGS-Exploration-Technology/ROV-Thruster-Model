%Thrust_to_Throttle_Table = csvread('Thrust_to_Throttle_Table.csv');
%Table_Thrusts = Thrust_to_Throttle_Table(:,1);
%Table_Throttles = Thrust_to_Throttle_Table(:,2);

Thrust = 19.5; %[N]

Throttle = thrust_to_throttle(Thrust, Table_Thrusts, Table_Throttles);

%Set up model params
n_samples = 10000;
Throttle_vector = Throttle*ones(1,n_samples);
Va = zeros(size(Throttle_vector));
dt = 0.001; %[s]
rho = 1027; %[kg/m^3] Density of seawater


%Thruster_Params = setup_thrusters_standard_ROV_full_vector_rev6;
%Thruster_Config = Thruster_Params(1); %Parameters for one thruster

load('Thruster_Config2.mat'); %Load previously fit thruster model params
Thruster_Config.RH_prop = 0;
Throttle_to_n_Table = csvread('Throttle_to_n_Table.csv');
Thruster_Config.Table_Throttles = Throttle_to_n_Table(:,1);
Thruster_Config.Table_n = Throttle_to_n_Table(:,2);

[t_fit, T_fit, Q_fit, n_fit] = sim_thruster_throttle(Throttle_vector, rho, Thruster_Config, dt);

%Plot results
figure;
subplot(3,1,1);
plot(t_fit, T_fit, '--k'); 
ylabel('Force [N]');
%axis([0 3 0 100]);
grid on;

subplot(3,1,2);
plot(t_fit, Q_fit, '--k');
ylabel('Torque [Nm]');
%axis([0 3 -1 20]);
grid on;

subplot(3,1,3);
plot(t_fit, n_fit, '--k');
ylabel('Prop Speed [rpm]');
grid on;
%axis([0 3 -1 35]);
