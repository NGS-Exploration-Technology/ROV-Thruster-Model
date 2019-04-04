function loss = fitness_fcn_rpm_thrust_torque(x,n_data,T_data,Q_data)
% Initialize optimization parameters:
cTn = x(1);
cQn = x(2);

% Calculate Thrust/Torque:
T = cTn*n_data.*abs(n_data); % [N] thrust
Q = cQn*n_data.*abs(n_data); % [Nm] torque

% Calculate cost function:
loss = sqrt(mean((T-T_data).^2)+mean((Q-Q_data).^2));
end