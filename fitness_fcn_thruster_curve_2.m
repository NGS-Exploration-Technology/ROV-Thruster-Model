function loss = fitness_fcn_thruster_curve_2(x,t,tv,n,v,u,constants)
% Fitness function for 2nd order rpm/velocity model in water

% Preallocate constants:
kv1 = constants{1};
kv2 = constants{2};
dv1 = constants{3};
dv2 = constants{4};
kV = constants{5};
kVV = constants{6};
kTV = constants{7};
cTn = constants{8};
cTnv = constants{9};
cQn = constants{10};
cQnv = constants{11};
dt = constants{12};
dtv = constants{13};
kn = x(1);
kq = x(2);
kv = x(3);
kvv = x(4);
kt = x(5);

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
tvp3 = tv{1};
tvp2 = tv{2};
tvm2 = tv{3};
tvm3 = tv{4};
vp3 = v{1};
vp2 = v{2};
vm2 = v{3};
vm3 = v{4};

% Calculate model:
n_fitp3 = zeros(length(tp3),1);
n_fitp2 = zeros(length(tp2),1);
n_fitm2 = zeros(length(tm2),1);
n_fitm3 = zeros(length(tm3),1);
v_fitp3 = zeros(length(tvp3),1);
v_fitp2 = zeros(length(tvp2),1);
v_fitm2 = zeros(length(tvm2),1);
v_fitm3 = zeros(length(tvm3),1);
vap3 = zeros(length(tvp3),1);
vap2 = zeros(length(tvp2),1);
vam2 = zeros(length(tvm2),1);
vam3 = zeros(length(tvm3),1);
idx = 2;
for i = 2:length(n_fitp3)
    Q = cQn*n_fitp3(i-1)*abs(n_fitp3(i-1))-cQnv*v_fitp3(idx-1)*abs(n_fitp3(i-1));
    d_n = -kn*n_fitp3(i-1)-kq*Q+kv2*(up3(i-1)-dv2);
    n_fitp3(i) = n_fitp3(i-1)+d_n*dt;
    
    if tvp3(idx-1) <= tp3(i-1)
        T = cTn*n_fitp3(i-1)*abs(n_fitp3(i-1))-cTnv*v_fitp3(idx-1)*abs(n_fitp3(i-1));
        d_v = -kv*v_fitp3(idx-1)-kvv*(v_fitp3(idx-1)-vap3(idx-1))*abs(v_fitp3(idx-1))+kt*T;
        d_va = -kV*vap3(idx-1)-kVV*vap3(idx-1)*abs(vap3(idx-1))+kTV*T;
        v_fitp3(idx) = v_fitp3(idx-1)+d_v*dtv;
        vap3(idx) = vap3(idx-1)+d_va*dtv;
        idx = idx+1;
    end
end
idx = 2;
for i = 2:length(n_fitp2)
    Q = cQn*n_fitp2(i-1)*abs(n_fitp2(i-1))-cQnv*v_fitp2(idx-1)*abs(n_fitp2(i-1));
    d_n = -kn*n_fitp2(i-1)-kq*Q+kv2*(up2(i-1)-dv2);
    n_fitp2(i) = n_fitp2(i-1)+d_n*dt;
    
    if tvp2(idx-1) <= tp2(i-1)
        T = cTn*n_fitp2(i-1)*abs(n_fitp2(i-1))-cTnv*v_fitp2(idx-1)*abs(n_fitp2(i-1));
        d_v = -kv*v_fitp2(idx-1)-kvv*(v_fitp2(idx-1)-vap2(idx-1))*abs(v_fitp2(idx-1))+kt*T;
        d_va = -kV*vap2(idx-1)-kVV*vap2(idx-1)*abs(vap2(idx-1))+kTV*T;
        v_fitp2(idx) = v_fitp2(idx-1)+d_v*dtv;
        vap2(idx) = vap2(idx-1)+d_va*dtv;
        idx = idx+1;
    end
end
idx = 2;
for i = 2:length(n_fitm2)
    Q = cQn*n_fitm2(i-1)*abs(n_fitm2(i-1))-cQnv*v_fitm2(idx-1)*abs(n_fitm2(i-1));
    d_n = -kn*n_fitm2(i-1)-kq*Q+kv1*(um2(i-1)-dv1);
    n_fitm2(i) = n_fitm2(i-1)+d_n*dt;
    
    if tvm2(idx-1) <= tm2(i-1)
        T = cTn*n_fitm2(i-1)*abs(n_fitm2(i-1))-cTnv*v_fitm2(idx-1)*abs(n_fitm2(i-1));
        d_v = -kv*v_fitm2(idx-1)-kvv*(v_fitm2(idx-1)-vam2(idx-1))*abs(v_fitm2(idx-1))+kt*T;
        d_va = -kV*vam2(idx-1)-kVV*vam2(idx-1)*abs(vam2(idx-1))+kTV*T;
        v_fitm2(idx) = v_fitm2(idx-1)+d_v*dtv;
        vam2(idx) = vam2(idx-1)+d_va*dtv;
        idx = idx+1;
    end
end
idx = 2;
for i = 2:length(n_fitm3)
    Q = cQn*n_fitm3(i-1)*abs(n_fitm3(i-1))-cQnv*v_fitm3(idx-1)*abs(n_fitm3(i-1));
    d_n = -kn*n_fitm3(i-1)-kq*Q+kv1*(um3(i-1)-dv1);
    n_fitm3(i) = n_fitm3(i-1)+d_n*dt;
    
    if tvm3(idx-1) <= tm3(i-1)
        T = cTn*n_fitm3(i-1)*abs(n_fitm3(i-1))-cTnv*v_fitm3(idx-1)*abs(n_fitm3(i-1));
        d_v = -kv*v_fitm3(idx-1)-kvv*(v_fitm3(idx-1)-vam3(idx-1))*abs(v_fitm3(idx-1))+kt*T;
        d_va = -kV*vam3(idx-1)-kVV*vam3(idx-1)*abs(vam3(idx-1))+kTV*T;
        v_fitm3(idx) = v_fitm3(idx-1)+d_v*dtv;
        vam3(idx) = vam3(idx-1)+d_va*dtv;
        idx = idx+1;
    end
end

loss = sqrt(mean((np3-n_fitp3).^2))+sqrt(mean((np2-n_fitp2).^2))+sqrt(mean((nm2-n_fitm2).^2))+sqrt(mean((nm3-n_fitm3).^2))+sqrt(mean((vm2(1:110)-v_fitm2(1:110)).^2))+sqrt(mean((vm2(179:end)-v_fitm2(179:end)).^2))+sqrt(mean((vm3(1:68)-v_fitm3(1:68)).^2))+sqrt(mean((vm3(125:end)-v_fitm3(125:end)).^2));
end