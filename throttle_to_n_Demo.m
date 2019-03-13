Throttle = -5; %[V]

Throttle_to_n_Table = csvread('Throttle_to_n_Table.csv');
Table_Throttles = Throttle_to_n_Table(:,1);
Table_n = Throttle_to_n_Table(:,2);
n=throttle_to_n(Throttle, Table_Throttles, Table_n) %[rpm]
n_Interp = interp1(Table_Throttles, Table_n, Throttle, 'linear') %[rpm]
