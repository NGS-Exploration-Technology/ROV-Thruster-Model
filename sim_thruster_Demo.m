Thrust_to_Throttle_Table = csvread('Thrust_to_Throttle_Table.csv');
Table_Thrusts = Thrust_to_Throttle_Table(:,1);
Table_Throttles = Thrust_to_Throttle_Table(:,2);

%Thrust = 0; %[N]

%Throttle = 5*thrust_to_throttle(Thrust, Table_Thrusts, Table_Throttles);

Throttle = -1;

%Set up model params
n_samples = 10000;
Throttle_vector = Throttle*ones(1,n_samples);
Va = zeros(size(Throttle_vector));
dt = 0.001; %[s]
rho = 1027; %[kg/m^3] Density of seawater

load('Thruster_Config.mat'); %Load previously fit thruster model params

[t_fit, T_fit, Q_fit, n_fit, u_fit] = sim_thruster(Va, Throttle_vector , rho, Thruster_Config, dt);

%Plot results
figure;
subplot(4,1,1);
plot(t_fit, T_fit, '--k'); 
ylabel('Force [N]');
legend('Experimental', 'Model');
%axis([0 3 0 100]);
grid on;

subplot(4,1,2);
plot(t_fit, Q_fit, '--k');
ylabel('Torque [Nm]');
%axis([0 3 -1 20]);
grid on;
subplot(4,1,3);
plot(t_fit, n_fit, '--k');
ylabel('Prop Speed [rpm]');
%axis([0 3 -1 35]);
grid on;
subplot(4,1,4);
plot(t_fit,u_fit,'--k')
xlabel('Time [s]'); ylabel('Flow Speed [m/s]');
grid on;