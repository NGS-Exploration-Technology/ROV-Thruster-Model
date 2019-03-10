Thrust_to_Throttle_Table = csvread('Thrust_to_Throttle_Table.csv');
Table_Thrusts = Thrust_to_Throttle_Table(:,1);
Table_Throttles = Thrust_to_Throttle_Table(:,2);

%Thrust = 0; %[N]

%Throttle = 5*thrust_to_throttle(Thrust, Table_Thrusts, Table_Throttles);

n_command = -2305;

%Set up model params
n_samples = 10000;
n_command_vector = n_command*ones(1,n_samples);
Va = zeros(size(n_command_vector));
dt = 0.001; %[s]
rho = 1027; %[kg/m^3] Density of seawater

load('Thruster_Config2.mat'); %Load previously fit thruster model params
Thruster_Config.RH_prop = 0;

[t_fit, T_fit, Q_fit, n_fit] = sim_thruster(n_command_vector , rho, Thruster_Config, dt);

%Plot results
figure;
subplot(3,1,1);
plot(t_fit, T_fit, '--k'); 
ylabel('Force [N]');
legend('Experimental', 'Model');
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
%axis([0 3 -1 35]);
