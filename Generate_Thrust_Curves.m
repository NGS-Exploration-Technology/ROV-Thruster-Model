function x = Generate_Thrust_Curves()

%Set fminsearch options
options.Display= 'final';
options.MaxFunEvals = [];
options.MaxIter = [];
options.TolFun = 1E-12;
options.TolX = 1E-12;
options.FunValCheck= [];
options.OutputFcn = [];
options.PlotFcns = {@optimplotfval,@optimplotx};

% Read in data:
A = csvread('Tach_RPM_Thrust_Data.csv');
n = A(:,7);
T = A(:,5);
Q = A(:,6);

% fminsearch operation:
f = @(x)fitness_fcn_rpm_thrust_torque(x,n,T,Q);

[x,fval,exitflag,opt_output] = fminsearch(f,[1 1 -1 1 1 1 -1 1], options)

cT1 = x(1);
cT2 = x(2);
dT1 = x(3);
dT2 = x(4);
cQ1 = x(5);
cQ2 = x(6);
dQ1 = x(7);
dQ2 = x(8);

% Initialize constants:
D = 0.1151; %[m] propellor diameter
rho = 1027; %[kg/m^3] Density of seawater

% Precalculate quadratic term:
n2 = n.*abs(n);

% Delete problematic data:
n2Q = n2;
n2Q(7) = [];
Q(7) = [];

% Calculate deadband piecewise quadratic function:
alpha1 = [cT1*(n2(n2<=dT1)-dT1);0*n2(n2<dT2 & n2>dT1);cT2*(n2(n2>=dT2)-dT2)];
beta1 = [cQ1*(n2Q(n2Q<=dQ1)-dQ1);0*n2Q(n2Q<dQ2 & n2Q>dQ1);cQ2*(n2Q(n2Q>=dQ2)-dQ2)];

% Calculate Thrust/Torque:
T_fit = rho*D^4*alpha1; % [N] thrust
Q_fit = rho*D^5*beta1; % [Nm] torque

% figure
% subplot(2,1,1)
% plot(n,T,'o',n,T_fit,'--k')
% ylabel('Thrust [N]')
% legend('Experimental','Model','Location','Northwest')
% grid
% 
% n(7) = [];
% subplot(2,1,2)
% plot(n,-Q,'o',n,-Q_fit,'--k')
% xlabel('Propeller Speed [rpm]')
% ylabel('Torque [Nm]')
% grid

font = 12;
width = 1.5;
figure
yyaxis left
plot(n,T,'o',n,T_fit,'--b','LineWidth',width)
xlabel('Propeller Speed [rpm]','FontSize',font,'FontName','Times New Roman')
ylabel('Thrust [N]','FontSize',font,'FontName','Times New Roman')

n(7) = [];
yyaxis right
plot(n,-Q,'o',n,-Q_fit,'--r','LineWidth',width)
ylabel('Torque [Nm]','FontSize',font,'FontName','Times New Roman')

legend({'Experiment','Model','Experiment','Model'},'FontSize',font,'FontName','Times New Roman','Location','Northwest')
set(gca,'FontSize',font,'FontName','Times New Roman')
grid
end