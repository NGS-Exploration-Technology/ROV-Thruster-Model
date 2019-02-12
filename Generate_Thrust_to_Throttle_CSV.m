[x, Throttle_Voltage, n, T, Q] = Generate_Thrust_Curves();

%Fix non-unique zero points in thrust
T(9) = -1e-9;
T(11) = 1e-9;

%Prune uneeded zero points in thrust and prune corresponding vectors
T(10) = [];
Throttle_Voltage(10) = [];
Q(10) = [];
n(10) = [];

%Rescale Throttle Voltage to a -1 to 1 Throttle signal
Throttle = Throttle_Voltage./5;

%Throttles_Table = linspace(-1,1,128)';
%Thrusts_Table = interp1(Throttle,T,Throttles_Table,'makima');
Thrusts_Table = linspace(min(T),max(T),128)';
Throttles_Table = interp1(T,Throttle,Thrusts_Table,'pchip');
Torques_Table = interp1(T,Q,Thrusts_Table,'makima');

subplot(2,1,1);
plot(T,Throttle,'xb',Thrusts_Table, Throttles_Table,'.k');
subplot(2,1,2);
plot(T,Q,'xb',Thrusts_Table,Torques_Table,'.k');

filename = 'Thrust_to_Throttle_Table.csv';
Thrust_to_Throttle_Table = [Thrusts_Table Throttles_Table];
csvwrite(filename,Thrust_to_Throttle_Table);

filename = 'Thrust_to_Torque_Table.csv';
Thrust_to_Torque_Table = [Thrusts_Table Torques_Table];
csvwrite(filename, Thrust_to_Torque_Table);
