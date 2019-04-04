function loss = fitness_fcn_rpm_voltage(x,t,n,v,dt)
% Fitness function for 1st order rpm model in air

% Preallocate constants:
kn = abs(x(1));
kv2 = abs(x(2));
dv2 = .9254;

% Calculate model:
n_fit = zeros(length(t),1);
for i = 2:length(n_fit)
    d_n = -kn*n_fit(i-1)+kv2*(v(i-1)-dv2);
    n_fit(i) = n_fit(i-1)+d_n*dt;
end

loss = sqrt(mean((n-n_fit).^2));
end