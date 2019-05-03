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
[cTn,cQn] = Generate_Thrust_Curves();
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
Datap2 = csvread('2statep2V.csv');
tp2 = Datap2(1:length(lenp2),2);
np2 = 49.9723*Datap2(lenp2,21);
Tp2 = -Datap2(lenp2,16);
Qp2 = Datap2(lenp2,19);
up2 = [(linspace(dv2,2,1200)).';2*ones(length(lenp2)-1200,1)];
Datam2 = csvread('2statem2V.csv');
tm2 = Datam2(1:length(lenm2),2);
nm2 = -49.9723*Datam2(lenm2,21);
Tm2 = -Datam2(lenm2,16)-1.8;
Qm2 = Datam2(lenm2,19);
um2 = [(linspace(dv1,-2,1200)).';-2*ones(length(lenm2)-1200,1)];
Datam3 = csvread('2statem3V.csv');
tm3 = Datam3(1:length(lenm3),2);
nm3 = -49.9723*Datam3(lenm3,21);
Tm3 = -Datam3(lenm3,16);
Qm3 = Datam3(lenm3,19);
um3 = [(linspace(dv1,-3,1500)).';-3*ones(length(lenm3)-1500,1)];
dt = tm3(2)-tm3(1);

% figure
% plot(tm3,nm3)
% grid
% figure
% plot(tm3,Tm3)
% grid
% figure
% plot(tm3,Qm3)
% grid

t = {tp3,tp2,tm2,tm3};
n = {np3,np2,nm2,nm3};
u = {up3,up2,um2,um3};

%define Fitness Function
f = @(x)fitness_fcn_thruster_curve(x,t,n,u,dt,cQn,kv1,kv2,dv1,dv2);

[x fval exitflag opt_output] = fminsearch(f,[kn .0038], options)

kn = x(1);
kq = x(2);

% Save results:
Thruster_Config.kv1 = kv1;
Thruster_Config.kv2 = kv2;
Thruster_Config.dv1 = dv1;
Thruster_Config.dv2 = dv2;
Thruster_Config.cTn = cTn;
Thruster_Config.cQn = cQn;
Thruster_Config.kn = kn;
Thruster_Config.kq = kq;

% kq = ((59.44/110.1)*kv2*(3-dv2)-kv2*(2-dv2))/(59.44*110-59.44^2);
% kn = (kv2*(3-dv2)-kq*110.1^2)/110.1;

% Simulate model:
n_fitp3 = zeros(length(tp3),1);
n_fitp2 = zeros(length(tp2),1);
n_fitm2 = zeros(length(tm2),1);
n_fitm3 = zeros(length(tm3),1);
for i = 2:length(n_fitp3)
    d_n = -kn*n_fitp3(i-1)-kq*cQn*n_fitp3(i-1)*abs(n_fitp3(i-1))+kv2*(up3(i-1)-dv2);
    n_fitp3(i) = n_fitp3(i-1)+d_n*dt;
end
T_fitp3 = cTn*n_fitp3.*abs(n_fitp3);
Q_fitp3 = cQn*n_fitp3.*abs(n_fitp3);
for i = 2:length(n_fitp2)
    d_n = -kn*n_fitp2(i-1)-kq*cQn*n_fitp2(i-1)*abs(n_fitp2(i-1))+kv2*(up2(i-1)-dv2);
    n_fitp2(i) = n_fitp2(i-1)+d_n*dt;
end
T_fitp2 = cTn*n_fitp2.*abs(n_fitp2);
Q_fitp2 = cQn*n_fitp2.*abs(n_fitp2);
for i = 2:length(n_fitm2)
    d_n = -kn*n_fitm2(i-1)-kq*cQn*n_fitm2(i-1)*abs(n_fitm2(i-1))+kv1*(um2(i-1)-dv1);
    n_fitm2(i) = n_fitm2(i-1)+d_n*dt;
end
T_fitm2 = cTn*n_fitm2.*abs(n_fitm2);
Q_fitm2 = cQn*n_fitm2.*abs(n_fitm2);
for i = 2:length(n_fitm3)
    d_n = -kn*n_fitm3(i-1)-kq*cQn*n_fitm3(i-1)*abs(n_fitm3(i-1))+kv1*(um3(i-1)-dv1);
    n_fitm3(i) = n_fitm3(i-1)+d_n*dt;
end
T_fitm3 = cTn*n_fitm3.*abs(n_fitm3);
Q_fitm3 = cQn*n_fitm3.*abs(n_fitm3);

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
plot(tp3,Tp3,tp3,T_fitp3,'--k')
grid
subplot(4,1,2)
plot(tp2,Tp2,tp2,T_fitp2,'--k')
legend('Data','Model')
ylabel('Thrust, T, [N]')
grid
subplot(4,1,3)
plot(tm2,Tm2,tm2,T_fitm2,'--k')
grid
subplot(4,1,4)
plot(tm3,Tm3,tm3,T_fitm3,'--k')
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