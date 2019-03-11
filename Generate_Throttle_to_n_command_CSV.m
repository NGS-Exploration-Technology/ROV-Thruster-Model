[x, Throttle_Voltage, n, T, Q] = Generate_Thrust_Curves();

%Fix non-unique zero points in thrust
T(9) = -1e-9;
T(11) = 1e-9;

%Prune uneeded zero points in thrust and prune corresponding vectors
T(10) = [];
Throttle_Voltage(10) = [];
Q(10) = [];
n(10) = [];

%Rescale Throttle if necessary
Throttle = Throttle_Voltage;

%Throttles_Table = linspace(-1,1,128)';
%Thrusts_Table = interp1(Throttle,T,Throttles_Table,'makima');
Throttles_Table = linspace(-5,5,128)';
n_Table = interp1(Throttle,n,Throttles_Table,'makima');

plot(Throttle, n,'xb',Throttles_Table, n_Table,'.k');
grid on;
xlabel('Throttle [V]');
ylabel('Propellor Angular Velocity [rpm]');

filename = 'Throttle_to_n_Table.csv';
Throttle_to_n_Table = [Throttles_Table n_Table];
csvwrite(filename,Throttle_to_n_Table);
