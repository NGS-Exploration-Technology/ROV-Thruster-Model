function loss = fitness_fcn_thruster_curve(x,t,n,u,dt,cQn,kv1,kv2,dv1,dv2)
% Fitness function for 1st order rpm model in water

% Preallocate constants:
kn = x(1);
kq = x(2);

% Extract data:
tp3 = t{1};
tp2 = t{2};
tm2 = t{3};
tm3 = t{4};
np3 = n{1};
np2 = n{2};
nm2 = n{3};
nm3 = n{4};
up3 = u{1};
up2 = u{2};
um2 = u{3};
um3 = u{4};

% Calculate model:
n_fitp3 = zeros(length(tp3),1);
n_fitp2 = zeros(length(tp2),1);
n_fitm2 = zeros(length(tm2),1);
n_fitm3 = zeros(length(tm3),1);
for i = 2:length(n_fitp3)
    d_n = -kn*n_fitp3(i-1)-kq*cQn*n_fitp3(i-1)*abs(n_fitp3(i-1))+kv2*(up3(i-1)-dv2);
    n_fitp3(i) = n_fitp3(i-1)+d_n*dt;
end
for i = 2:length(n_fitp2)
    d_n = -kn*n_fitp2(i-1)-kq*cQn*n_fitp2(i-1)*abs(n_fitp2(i-1))+kv2*(up2(i-1)-dv2);
    n_fitp2(i) = n_fitp2(i-1)+d_n*dt;
end
for i = 2:length(n_fitm2)
    d_n = -kn*n_fitm2(i-1)-kq*cQn*n_fitm2(i-1)*abs(n_fitm2(i-1))+kv1*(um2(i-1)-dv1);
    n_fitm2(i) = n_fitm2(i-1)+d_n*dt;
end
for i = 2:length(n_fitm3)
    d_n = -kn*n_fitm3(i-1)-kq*cQn*n_fitm3(i-1)*abs(n_fitm3(i-1))+kv1*(um3(i-1)-dv1);
    n_fitm3(i) = n_fitm3(i-1)+d_n*dt;
end

loss = sqrt(mean((np3-n_fitp3).^2))+sqrt(mean((np2-n_fitp2).^2))+sqrt(mean((nm2-n_fitm2).^2))+sqrt(mean((nm3-n_fitm3).^2));
end