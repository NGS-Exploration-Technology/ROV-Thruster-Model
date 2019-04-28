function x = optimize_ambient_dynamics()
% Read in ambient velocity data and fit to first order model

% Set fminsearch options:
options.Display= 'final';
options.MaxFunEvals = [];
options.MaxIter = [];
options.TolFun = 1E-12;
options.TolX = 1E-12;
options.FunValCheck= [];
options.OutputFcn = [];
options.PlotFcns = {@optimplotfval,@optimplotx};

% Initialize some parameters:
lenm3 = 123:358;
lenm2 = 129:354;
lenp2 = 75:301;
lenp3 = 90:323;

% Read in data:
Vrawm3 = load('Ambientm3V.dat');
Vm3 = -Vrawm3(lenm3,4);
tm3 = Vrawm3(lenm3,1)-Vrawm3(lenm3(1),1);
Tm3 = -[6.6444*tm3(1:38).^2;13*ones(length(lenm3)-38,1)];
Vrawm2 = load('Ambientm2V.dat');
Vm2 = -Vrawm2(lenm2,4);
tm2 = Vrawm2(lenm2,1)-Vrawm2(lenm2(1),1);
Tm2 = -[3.5*tm2(1:30).^2;4*ones(length(lenm2)-30,1)];
Vrawp2 = load('Ambientp2V.dat');
Vp2 = Vrawp2(lenp2,4);
tp2 = Vrawp2(lenp2,1)-Vrawp2(lenp2(1),1);
Tp2 = [3*tp2(1:25).^2;3.25*ones(length(lenp2)-25,1)];
Vrawp3 = load('Ambientp3V.dat');
Vp3 = Vrawp3(lenp3,4);
tp3 = Vrawp3(lenp3,1)-Vrawp3(lenp3(1),1);
Tp3 = [7*tp3(1:37).^2;11.3*ones(length(lenp3)-37,1)];
dt = tp3(2)-tp3(1);

t = {tp3,tp2,tm2,tm3};
v = {Vp3,Vp2,Vm2,Vm3};
T = {Tp3,Tp2,Tm2,Tm3};

% Traw = csvread('Ambientp3V.csv');
% Tdata = Traw(3580:12880,16)+1.6;
% timeT = Traw(1:length(3580:12880),2);
% 
% figure
% plot(timeVp3,Vp3)
% xlabel('Time [sec]')
% ylabel('Ambient Velocity [m/s]')
% grid
% figure
% plot(timeT,Tdata,timeVp3,[7*timeVp3(1:37).^2;11.3*ones(length(lenp3)-37,1)])
% xlabel('Time [sec]')
% ylabel('Thrust [N]')
% grid
% figure
% plot(timeT,49.9723*Traw(3580:12880,21))
% xlabel('Time [sec]')
% ylabel('Propeller Velocity [rad/s]')
% grid

% fminsearch operation:
f = @(x)fitness_fcn_ambient(x,t,v,T,dt);

[x,fval,exitflag,opt_output] = fminsearch(f,[.2 33 .006], options)

% Retrieve constants:
kV = x(1);
kVV = x(2);
kTV = x(3);

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

figure
subplot(4,1,1)
plot(tp3,Vp3,tp3,v_fitp3,'--k')
legend('Data','Model')
grid
subplot(4,1,2)
plot(tp2,Vp2,tp2,v_fitp2,'--k')
ylabel('Ambient Velocity, V, [m/s]')
grid
subplot(4,1,3)
plot(tm2,Vm2,tm2,v_fitm2,'--k')
grid
subplot(4,1,4)
plot(tm3,Vm3,tm3,v_fitm3,'--k')
xlabel('Time, t, [sec]')
grid

end