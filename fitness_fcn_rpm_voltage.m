function loss = fitness_fcn_rpm_voltage(x,t,n,v,dt)
% Fitness function for 1st order rpm model in air

% Preallocate constants:
kn = x(1);
kv1 = x(2);
kv2 = x(3);
dv1 = -.83;
dv2 = .86;

% Extract data:
tp5 = t{1};
tp3 = t{2};
tm3 = t{3};
tm5 = t{4};
np5 = n{1};
np3 = n{2};
nm3 = n{3};
nm5 = n{4};
vp5 = v{1};
vp3 = v{2};
vm3 = v{3};
vm5 = v{4};

% Calculate model:
n_fitp5 = zeros(length(tp5),1);
n_fitp3 = zeros(length(tp3),1);
n_fitm3 = zeros(length(tm3),1);
n_fitm5 = zeros(length(tm5),1);
for i = 2:length(n_fitp5)
    d_n = -kn*n_fitp5(i-1)+kv2*(vp5(i-1)-dv2);
    n_fitp5(i) = n_fitp5(i-1)+d_n*dt;
end
for i = 2:length(n_fitp3)
    d_n = -kn*n_fitp3(i-1)+kv2*(vp3(i-1)-dv2);
    n_fitp3(i) = n_fitp3(i-1)+d_n*dt;
end
for i = 2:length(n_fitm3)
    d_n = -kn*n_fitm3(i-1)+kv1*(vm3(i-1)-dv1);
    n_fitm3(i) = n_fitm3(i-1)+d_n*dt;
end
for i = 2:length(n_fitm5)
    d_n = -kn*n_fitm5(i-1)+kv1*(vm5(i-1)-dv1);
    n_fitm5(i) = n_fitm5(i-1)+d_n*dt;
end

loss = sqrt(mean((np5-n_fitp5).^2))+sqrt(mean((np3-n_fitp3).^2))+sqrt(mean((nm3-n_fitm3).^2))+sqrt(mean((nm5-n_fitm5).^2));
end