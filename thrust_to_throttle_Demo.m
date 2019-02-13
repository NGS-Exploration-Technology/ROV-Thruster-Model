Thrust_to_Throttle_Table = csvread('Thrust_to_Throttle_Table.csv');
Thrust_to_Torque_Table = csvread('Thrust_to_Torque_Table.csv');
Table_Thrusts = Thrust_to_Throttle_Table(:,1);
Table_Throttles = Thrust_to_Throttle_Table(:,2);
Table_Torques = Thrust_to_Torque_Table(:,2);

Desired_Thrust = 0; %[N]

Throttle = thrust_to_throttle(Desired_Thrust, Table_Thrusts, Table_Throttles) %[0-1]
Throttle_Interp = interp1(Table_Thrusts, Table_Throttles, Desired_Thrust, 'pchip') %[0-1]

Torque = thrust_to_torque(Desired_Thrust, Table_Thrusts, Table_Torques) %[Nm]
Torque_Interp = interp1(Table_Thrusts, Table_Torques, Desired_Thrust, 'pchip') %[Nm]
