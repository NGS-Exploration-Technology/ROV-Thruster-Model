%% Average Thrust Histories
%% Initialize
clc; clear; close all;

%% Run Routine
lpFilt = designfilt('lowpassiir','FilterOrder',4, ...
         'PassbandFrequency',15,'PassbandRipple',0.01, ...
         'StopbandAttenuation', 120, 'SampleRate',1e3);

Set1 = csvread('Fwd_Rvrs_Set1.csv');
Set2 = csvread('Fwd_Rvrs_Set2.csv');
Set3 = csvread('Fwd_Rvrs_Set3.csv');
Set4 = csvread('Fwd_Rvrs_Set4.csv');
Set5 = csvread('Fwd_Rvrs_Set5.csv');

len1 = 770:5270;
len2 = 2440:6940;
len3 = 810:5310;
len4 = 200:4700;
len5 = 910:5410;

t1 = Set1(len1,2);
t2 = Set2(len2,2);
t3 = Set3(len3,2);
t4 = Set4(len4,2);
t5 = Set5(len5,2);

n1 = (477.43*Set1(len1,21)+5.3255)*2*pi/60;
n2 = (477.43*Set2(len2,21)+5.3255)*2*pi/60;
n3 = (477.43*Set3(len3,21)+5.3255)*2*pi/60;
n4 = (477.43*Set4(len4,21)+5.3255)*2*pi/60;
n5 = (477.43*Set5(len5,21)+5.3255)*2*pi/60;

F1 = Set1(len1,16)*-1;
F2 = Set2(len2,16)*-1;
F3 = Set3(len3,16)*-1;
F4 = Set4(len4,16)*-1;
F5 = Set5(len5,16)*-1;

% figure
% yyaxis left
% plot(t1,F1,t2,F2,t3,F3,t4,F4,t5,F5)
% axis([.5 5.5 -20 20])
% yyaxis right
% plot(t5,n5)
% grid

Fmat = [F1.';F2.';F3.';F4.';F5.'];
Favg = mean(Fmat);
nmat = [n1.';n2.';n3.';n4.';n5.'];
navg = mean(nmat);
Flpf = filter(lpFilt,Favg);

figure
plot(navg,Favg)
grid

% Setup simulation parameters
thruster_config;
rho = 1027;
Diam = Thruster_Config.D;
L = Thruster_Config.L;
kn = Thruster_Config.kn1;
kq = Thruster_Config.kq*2*pi/60;
kv = 1/.05;
kvv = 2/L;
kt = (rho*L*pi*(Diam/2)^2)^-1;
% kv1 = 7578*2*pi/60;
kv2 = 6770*2*pi/60;
kv1 = kv2;
dv1 = -.8475;
dv2 = .9254;
cTn = .004243*1.725*1.8*.5;
cTv = -9.435*30;
cQn = .000626*3;
cQv = -1.212*10;
udotmax = 1/.01;

% Customize ICs and simulation properties:
dt = .001;
t = 0:dt:4.5;
n = zeros(length(t),1);
v = zeros(length(t),1);
va = zeros(length(t),1);
u = [linspace(0,3,500).';3*ones(1500,1);linspace(3,-3,1000).';-3*ones(1000,1);linspace(-3,0,501).'];
T = zeros(length(t),1);
Q = zeros(length(t),1);
n(1) = 0;
v(1) = 0;

% Run propagation:
for i = 1:(length(t)-1)
    T(i) = rho*Diam^4*(cTn*n(i)*abs(n(i))+cTv*v(i)*abs(v(i)));
    Q(i) = rho*Diam^5*(cQn*n(i)*abs(n(i))+cQv*v(i)*abs(v(i)));

    if i == 1
        uprev = u(i);
    else
        uprev = u(i-1);
    end
    if (abs(u(i)-uprev)/dt)>udotmax
        u(i) = uprev+sign(u(i)-uprev)*udotmax*dt;
    end
    if u(i) <= dv1
        du = kv1*(u(i)-dv1);
    elseif u(i) >= dv2
        du = kv2*(u(i)-dv2);
    else
        du = 0;
    end

    ndot = -kn*n(i)-kq*Q(i)+du; % note sign on Q is opposite to what we've done before
    vdot = -kv*v(i)-kvv*abs(v(i))*(v(i)-va(i))+kt*T(i);
    n(i+1) = n(i)+ndot*dt;
    v(i+1) = v(i)+vdot*dt;
end
T(end) = rho*Diam^4*(cTn*n(end)*abs(n(end))+cTv*v(end)*abs(v(end)));
Q(end) = rho*Diam^5*(cQn*n(end)*abs(n(end))+cQv*v(end)*abs(v(end)));

n_old = -149:.1:150;
hold on
plot(n,T,n,rho*Diam^4*(cTn*n.*abs(n)-1.48693*.7*abs(n).*v),n_old,rho*Diam^4*.7*Thruster_Config.cT1/((2*pi/60)^2)*n_old.*abs(n_old),'--k')
xlabel('Propeller Velocity, n, [rad/s]')
ylabel('Forward Thrust, T, [N]')
legend('Experiment','Two-State','Fossen','Single-State','Location','Northwest')
hold off