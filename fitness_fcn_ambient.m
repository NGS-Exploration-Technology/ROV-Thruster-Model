function loss = fitness_fcn_ambient(x,t,v,T,dt)
% Fitness function for 1st order ambient velocity model

% Preallocate constants:
kV = x(1);
kVV = x(2);
kTV = x(3);

% Extract data:
tp3 = t{1};
tp2 = t{2};
tm2 = t{3};
tm3 = t{4};
vp3 = v{1};
vp2 = v{2};
vm2 = v{3};
vm3 = v{4};
Tp3 = T{1};
Tp2 = T{2};
Tm2 = T{3};
Tm3 = T{4};

% Calculate model:
v_fitp3 = zeros(length(tp3),1);
v_fitp2 = zeros(length(tp2),1);
v_fitm2 = zeros(length(tm2),1);
v_fitm3 = zeros(length(tm3),1);
for i = 2:length(v_fitp3)
    d_v = -kV*v_fitp3(i-1)-kVV*v_fitp3(i-1)*abs(v_fitp3(i-1))+kTV*Tp3(i-1);
    v_fitp3(i) = v_fitp3(i-1)+d_v*dt;
end
for i = 2:length(v_fitp2)
    d_v = -kV*v_fitp2(i-1)-kVV*v_fitp2(i-1)*abs(v_fitp2(i-1))+kTV*Tp2(i-1);
    v_fitp2(i) = v_fitp2(i-1)+d_v*dt;
end
for i = 2:length(v_fitm2)
    d_v = -kV*v_fitm2(i-1)-kVV*v_fitm2(i-1)*abs(v_fitm2(i-1))+kTV*Tm2(i-1);
    v_fitm2(i) = v_fitm2(i-1)+d_v*dt;
end
for i = 2:length(v_fitm3)
    d_v = -kV*v_fitm3(i-1)-kVV*v_fitm3(i-1)*abs(v_fitm3(i-1))+kTV*Tm3(i-1);
    v_fitm3(i) = v_fitm3(i-1)+d_v*dt;
end

loss = sqrt(mean((vp3-v_fitp3).^2))+sqrt(mean((vp2-v_fitp2).^2))+sqrt(mean((vm2-v_fitm2).^2))+sqrt(mean((vm3-v_fitm3).^2));
end