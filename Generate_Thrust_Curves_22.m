function [cTn,cTnv,cQn,cQnv] = Generate_Thrust_Curves_22()

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
A = csvread('Steady_State_Data.csv');
v = A(:,5);
n = A(:,4);
T = A(:,2);
Q = A(:,3);

% fminsearch operation:
f = @(x)sqrt(mean((x(1)*n.*abs(n)-x(2)*v.*abs(n)-x(3)*v.*abs(v)-T).^2)+mean((x(4)*n.*abs(n)-x(5)*v.*abs(n)-x(6)*v.*abs(v)-Q).^2));

[x,fval,exitflag,opt_output] = fminsearch(f,[.0024 1.0139 150.8548 2.8e-6 .00000000001 .0000001], options)

cTn = x(1);
cTnv = x(2);
cTv = x(3);
cQn = x(4);
cQnv = x(5);
cQv = x(6);

% Calculate Thrust/Torque:
T_fit = cTn*n.*abs(n)-cTnv*v.*abs(n)-cTv*v.*abs(v); % [N] thrust
Q_fit = cQn*n.*abs(n)-cQnv*v.*abs(n)-cQv*v.*abs(v); % [Nm] torque

% Generate Plots:
% rho = 1027;
% D = .1151;
% 
% font = 12;
% width = 1.5;
% figure
% yyaxis left
% plot(v./(n*D),T./(rho*D^4*n.*abs(n)),'o',v./(n*D),T_fit./(rho*D^4*n.*abs(n)),'--b','LineWidth',width)
% xlabel('Propeller Speed [rad/s]','FontSize',font,'FontName','Times New Roman')
% ylabel('Thrust [N]','FontSize',font,'FontName','Times New Roman')
% 
% yyaxis right
% plot(v./(n*D),Q./(rho*D^5*n.*abs(n)),'o',v./(n*D),Q_fit./(rho*D^5*n.*abs(n)),'--r','LineWidth',width)
% ylabel('Torque [Nm]','FontSize',font,'FontName','Times New Roman')
% 
% legend({'Experiment','Model','Experiment','Model'},'FontSize',font,'FontName','Times New Roman','Location','Northwest')
% set(gca,'FontSize',font,'FontName','Times New Roman')
% grid

font = 12;
width = 1.5;
figure
yyaxis left
plot(n,T,'o',n,T_fit,'--b','LineWidth',width)
xlabel('Propeller Speed [rad/s]','FontSize',font,'FontName','Times New Roman')
ylabel('Thrust [N]','FontSize',font,'FontName','Times New Roman')

yyaxis right
plot(n,Q,'o',n,Q_fit,'--r','LineWidth',width)
ylabel('Torque [Nm]','FontSize',font,'FontName','Times New Roman')

legend({'Experiment','Model','Experiment','Model'},'FontSize',font,'FontName','Times New Roman','Location','Northwest')
set(gca,'FontSize',font,'FontName','Times New Roman')
grid
end