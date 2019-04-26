function x = optimize_voltage_constant()
% Read in rpm data and fit to first order model

%Set fminsearch options
options.Display= 'final';
options.MaxFunEvals = [];
options.MaxIter = [];
options.TolFun = 1E-12;
options.TolX = 1E-12;
options.FunValCheck= [];
options.OutputFcn = [];
options.PlotFcns = {@optimplotfval,@optimplotx};

% Initialize some parameters:
lenp5 = 4635:6000;
lenp3 = 3185:5200;
lenm3 = 9714:11700;
lenm5 = 6800:8800;
dv1 = -.83;
dv2 = .86;

% Read in data:
Ap5 = csvread('201807201712_Air_RPM_Data.csv');
tp5 = Ap5(1:length(lenp5),2);
np5 = (477.43*Ap5(lenp5,21)+5.3255)*2*pi/60;
vp5 = 5*ones(length(lenp5),1);
Ap3 = csvread('201808061450_3V.csv');
tp3 = Ap3(1:length(lenp3),2);
np3 = (477.43*Ap3(lenp3,21)+5.3255)*2*pi/60;
vp3 = 3*ones(length(lenp3),1);
Am3 = csvread('201808061453_-3V.csv');
tm3 = Am3(1:length(lenm3),2);
nm3 = -(477.43*Am3(lenm3,21)+5.3255)*2*pi/60;
vm3 = -3*ones(length(lenm3),1);
Am5 = csvread('201808061456_-5V.csv');
tm5 = Am5(1:length(lenm5),2);
nm5 = -(477.43*Am5(lenm5,21)+5.3255)*2*pi/60;
vm5 = -5*ones(length(lenm5),1);
dt = tp5(2)-tp5(1);

np3mf = median_filter(np3,20);
nm3mf = median_filter(nm3,20);

t = {tp5,tp3,tm3,tm5};
n = {np5,np3,nm3,nm5};
v = {vp5,vp3,vm3,vm5};

plot(tm3,nm3,tm3,nm3mf)
grid

% fminsearch operation:
f = @(x)fitness_fcn_rpm_voltage(x,t,n,v,dt);

[x,fval,exitflag,opt_output] = fminsearch(f,[10 700 700], options)

% Set constants:
kn = x(1);
kv1 = x(2);
kv2 = x(3);

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

figure
plot(tp5,np5,tp5,n_fitp5,'--k')
legend('Data','Model')
xlabel('Time, t, [sec]')
ylabel('Rotor Speed, n, [rpm]')
grid
figure
plot(tp3,np3,tp3,np3mf,tp3,n_fitp3,'--k')
legend('Data','Filter','Model')
xlabel('Time, t, [sec]')
ylabel('Rotor Speed, n, [rpm]')
grid
figure
plot(tm3,nm3,tm3,nm3mf,tm3,n_fitm3,'--k')
legend('Data','Filter','Model')
xlabel('Time, t, [sec]')
ylabel('Rotor Speed, n, [rpm]')
grid
figure
plot(tm5,nm5,tm5,n_fitm5,'--k')
legend('Data','Model')
xlabel('Time, t, [sec]')
ylabel('Rotor Speed, n, [rpm]')
grid

end