clear; close all; clc;

%Set fminsearch options
options.Display= 'final';
options.MaxFunEvals = 5000;
options.MaxIter = [];
options.TolFun = 1E-12;
options.TolX = 1E-12;
options.FunValCheck= [];
options.OutputFcn = [];
options.PlotFcns = {@optimplotfval,@optimplotx};

% Define some parameters
[kn,kv1,kv2] = optimize_voltage_constant();
[kV,kVV,kTV] = optimize_ambient_dynamics();
[cTn,cTnv,cQn,cQnv] = Generate_Thrust_Curves_2();
lenp3 = 4300:13500;
lenp2 = 3450:12500;
lenm2 = 4550:13500;
lenm3 = 4500:14000;
dv1 = -.83;
dv2 = .86;

% Extract data
Datap3 = csvread('2statep3V.csv');
tp3 = Datap3(1:length(lenp3),2);
np3 = 49.9723*Datap3(lenp3,21);
Tp3 = -Datap3(lenp3,16);
Qp3 = Datap3(lenp3,19);
up3 = [(linspace(dv2,3,1600)).';3*ones(length(lenp3)-1600,1)];
Vrawp3 = load('2statep3V.dat');
vp3 = Vrawp3(109:339,4);
tvp3 = Vrawp3(109:339,1)-Vrawp3(109,1);
Datap2 = csvread('2statep2V.csv');
tp2 = Datap2(1:length(lenp2),2);
np2 = 49.9723*Datap2(lenp2,21);
Tp2 = -Datap2(lenp2,16);
Qp2 = Datap2(lenp2,19);
up2 = [(linspace(dv2,2,1200)).';2*ones(length(lenp2)-1200,1)];
Vrawp2 = load('2statep2V.dat');
vp2 = Vrawp2(87:314,4);
tvp2 = Vrawp2(87:314,1)-Vrawp2(87,1);
Datam2 = csvread('2statem2V.csv');
tm2 = Datam2(1:length(lenm2),2);
nm2 = -49.9723*Datam2(lenm2,21);
Tm2 = -Datam2(lenm2,16)-1.8;
Qm2 = Datam2(lenm2,19);
um2 = [(linspace(dv1,-2,1200)).';-2*ones(length(lenm2)-1200,1)];
Vrawm2 = load('2statem2V.dat');
vm2 = -Vrawm2(115:339,4);
tvm2 = Vrawm2(115:339,1)-Vrawm2(115,1);
Datam3 = csvread('2statem3V.csv');
tm3 = Datam3(1:length(lenm3),2);
nm3 = -49.9723*Datam3(lenm3,21);
Tm3 = -Datam3(lenm3,16);
Qm3 = Datam3(lenm3,19);
um3 = [(linspace(dv1,-3,1500)).';-3*ones(length(lenm3)-1500,1)];
Vrawm3 = load('2statem3V.dat');
vm3 = -Vrawm3(114:352,4);
tvm3 = Vrawm3(114:352,1)-Vrawm3(114,1);
dt = tm3(2)-tm3(1);
dtv = tvm3(2)-tvm3(1);

% figure
% plot(tm3,nm3)
% grid
% figure
% plot(tm3,Tm3)
% grid
% figure
% plot(tm3,Qm3)
% grid
% figure
% plot(time_vm3,vm3)
% grid

cTn = .55*cTn;
cTnv = .3*cTnv;

t = {tp3,tp2,tm2,tm3};
n = {np3,np2,nm2,nm3};
u = {up3,up2,um2,um3};
tv = {tvp3,tvp2,tvm2,tvm3};
v = {vp3,vp2,vm2,vm3};
constants = {kv1,kv2,dv1,dv2,kV,kVV,kTV,cTn,cTnv,cQn,cQnv,dt,dtv};

%define Fitness Function
f = @(x)fitness_fcn_thruster_curve_2(x,t,tv,n,v,u,constants);

[x fval exitflag opt_output] = fminsearch(f,[kn .0038 kV 10 kTV], options)

kn = x(1);
kq = x(2);
kv = x(3);
kvv = abs(x(4));
kt = x(5);

% Save results:
Thruster_Config.kv1 = kv1;
Thruster_Config.kv2 = kv2;
Thruster_Config.dv1 = dv1;
Thruster_Config.dv2 = dv2;
Thruster_Config.cTn = cTn;
Thruster_Config.cTnv = cTnv;
Thruster_Config.cQn = cQn;
Thruster_Config.cQnv = cQnv;
Thruster_Config.kn = kn;
Thruster_Config.kq = kq;
Thruster_Config.kv = kv;
Thruster_Config.kvv = kvv;
Thruster_Config.kt = kt;
Thruster_Config.kV = kV;
Thruster_Config.kVV = kVV;
Thruster_Config.kTV = kTV;

% Simulate model:
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
T_fitp3 = zeros(length(tvp3),1);
T_fitp2 = zeros(length(tvp2),1);
T_fitm2 = zeros(length(tvm2),1);
T_fitm3 = zeros(length(tvm3),1);
Q_fitp3 = zeros(length(tp3),1);
Q_fitp2 = zeros(length(tp2),1);
Q_fitm2 = zeros(length(tm2),1);
Q_fitm3 = zeros(length(tm3),1);
idx = 2;
for i = 2:length(n_fitp3)
    Q_fitp3(i-1) = cQn*n_fitp3(i-1)*abs(n_fitp3(i-1))-cQnv*v_fitp3(idx-1)*abs(n_fitp3(i-1));
    d_n = -kn*n_fitp3(i-1)-kq*Q_fitp3(i-1)+kv2*(up3(i-1)-dv2);
    n_fitp3(i) = n_fitp3(i-1)+d_n*dt;
    
    if tvp3(idx-1) <= tp3(i-1)
        T_fitp3(idx-1) = cTn*n_fitp3(i-1)*abs(n_fitp3(i-1))-cTnv*v_fitp3(idx-1)*abs(n_fitp3(i-1));
        d_v = -kv*v_fitp3(idx-1)-kvv*(v_fitp3(idx-1)-vap3(idx-1))*abs(v_fitp3(idx-1))+kt*T_fitp3(idx-1);
        d_va = -kV*vap3(idx-1)-kVV*vap3(idx-1)*abs(vap3(idx-1))+kTV*T_fitp3(idx-1);
        v_fitp3(idx) = v_fitp3(idx-1)+d_v*dtv;
        vap3(idx) = vap3(idx-1)+d_va*dtv;
        idx = idx+1;
    end
end
idx = 2;
for i = 2:length(n_fitp2)
    Q_fitp2(i-1) = cQn*n_fitp2(i-1)*abs(n_fitp2(i-1))-cQnv*v_fitp2(idx-1)*abs(n_fitp2(i-1));
    d_n = -kn*n_fitp2(i-1)-kq*Q_fitp2(i-1)+kv2*(up2(i-1)-dv2);
    n_fitp2(i) = n_fitp2(i-1)+d_n*dt;
    
    if tvp2(idx-1) <= tp2(i-1)
        T_fitp2(idx-1) = cTn*n_fitp2(i-1)*abs(n_fitp2(i-1))-cTnv*v_fitp2(idx-1)*abs(n_fitp2(i-1));
        d_v = -kv*v_fitp2(idx-1)-kvv*(v_fitp2(idx-1)-vap2(idx-1))*abs(v_fitp2(idx-1))+kt*T_fitp2(idx-1);
        d_va = -kV*vap2(idx-1)-kVV*vap2(idx-1)*abs(vap2(idx-1))+kTV*T_fitp2(idx-1);
        v_fitp2(idx) = v_fitp2(idx-1)+d_v*dtv;
        vap2(idx) = vap2(idx-1)+d_va*dtv;
        idx = idx+1;
    end
end
idx = 2;
for i = 2:length(n_fitm2)
    Q_fitm2(i-1) = cQn*n_fitm2(i-1)*abs(n_fitm2(i-1))-cQnv*v_fitm2(idx-1)*abs(n_fitm2(i-1));
    d_n = -kn*n_fitm2(i-1)-kq*Q_fitm2(i-1)+kv1*(um2(i-1)-dv1);
    n_fitm2(i) = n_fitm2(i-1)+d_n*dt;
    
    if tvm2(idx-1) <= tm2(i-1)
        T_fitm2(idx-1) = cTn*n_fitm2(i-1)*abs(n_fitm2(i-1))-cTnv*v_fitm2(idx-1)*abs(n_fitm2(i-1));
        d_v = -kv*v_fitm2(idx-1)-kvv*(v_fitm2(idx-1)-vam2(idx-1))*abs(v_fitm2(idx-1))+kt*T_fitm2(idx-1);
        d_va = -kV*vam2(idx-1)-kVV*vam2(idx-1)*abs(vam2(idx-1))+kTV*T_fitm2(idx-1);
        v_fitm2(idx) = v_fitm2(idx-1)+d_v*dtv;
        vam2(idx) = vam2(idx-1)+d_va*dtv;
        idx = idx+1;
    end
end
idx = 2;
for i = 2:length(n_fitm3)
    Q_fitm3(i-1) = cQn*n_fitm3(i-1)*abs(n_fitm3(i-1))-cQnv*v_fitm3(idx-1)*abs(n_fitm3(i-1));
    d_n = -kn*n_fitm3(i-1)-kq*Q_fitm3(i-1)+kv1*(um3(i-1)-dv1);
    n_fitm3(i) = n_fitm3(i-1)+d_n*dt;
    
    if tvm3(idx-1) <= tm3(i-1)
        T_fitm3(idx-1) = cTn*n_fitm3(i-1)*abs(n_fitm3(i-1))-cTnv*v_fitm3(idx-1)*abs(n_fitm3(i-1));
        d_v = -kv*v_fitm3(idx-1)-kvv*(v_fitm3(idx-1)-vam3(idx-1))*abs(v_fitm3(idx-1))+kt*T_fitm3(idx-1);
        d_va = -kV*vam3(idx-1)-kVV*vam3(idx-1)*abs(vam3(idx-1))+kTV*T_fitm3(idx-1);
        v_fitm3(idx) = v_fitm3(idx-1)+d_v*dtv;
        vam3(idx) = vam3(idx-1)+d_va*dtv;
        idx = idx+1;
    end
end

% Generate plots
figure
subplot(4,1,1)
plot(tp3,np3,tp3,n_fitp3,'--k')
grid
subplot(4,1,2)
plot(tp2,np2,tp2,n_fitp2,'--k')
legend('Data','Model')
ylabel('Rotor Speed, n, [rad/s]')
grid
subplot(4,1,3)
plot(tm2,nm2,tm2,n_fitm2,'--k')
grid
subplot(4,1,4)
plot(tm3,nm3,tm3,n_fitm3,'--k')
xlabel('Time, t, [sec]')
grid

figure
subplot(4,1,1)
plot(tvp3,vp3,tvp3,v_fitp3,'--k')
grid
subplot(4,1,2)
plot(tvp2,vp2,tvp2,v_fitp2,'--k')
legend('Data','Model')
ylabel('Axial Velocity, v, [m/s]')
grid
subplot(4,1,3)
plot(tvm2,vm2,tvm2,v_fitm2,'--k')
grid
subplot(4,1,4)
plot(tvm3,vm3,tvm3,v_fitm3,'--k')
xlabel('Time, t, [sec]')
grid

figure
subplot(4,1,1)
plot(tp3,Tp3,tvp3,T_fitp3,'--k')
grid
subplot(4,1,2)
plot(tp2,Tp2,tvp2,T_fitp2,'--k')
legend('Data','Model')
ylabel('Thrust, T, [N]')
grid
subplot(4,1,3)
plot(tm2,Tm2,tvm2,T_fitm2,'--k')
grid
subplot(4,1,4)
plot(tm3,Tm3,tvm3,T_fitm3,'--k')
xlabel('Time, t, [sec]')
grid

figure
subplot(4,1,1)
plot(tp3,Qp3,tp3,Q_fitp3,'--k')
grid
subplot(4,1,2)
plot(tp2,Qp2,tp2,Q_fitp2,'--k')
legend('Data','Model')
ylabel('Torque, Q, [N*m]')
grid
subplot(4,1,3)
plot(tm2,Qm2,tm2,Q_fitm2,'--k')
grid
subplot(4,1,4)
plot(tm3,Qm3,tm3,Q_fitm3,'--k')
xlabel('Time, t, [sec]')
grid