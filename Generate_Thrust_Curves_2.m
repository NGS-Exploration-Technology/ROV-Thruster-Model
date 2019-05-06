function [cTn,cTnv,cQn,cQnv] = Generate_Thrust_Curves_2()

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
f = @(x)sqrt(mean((x(1)*n.*abs(n)-x(2)*v.*abs(n)-T).^2)+mean((x(3)*n.*abs(n)-x(4)*v.*abs(n)-Q).^2));

[x,fval,exitflag,opt_output] = fminsearch(f,[.0033 1.9347 2.8e-6 .00000000001], options)

cTn = x(1);
cTnv = x(2);
cQn = x(3);
cQnv = x(4);

% Calculate Thrust/Torque:
T_fit = cTn*n.*abs(n)-cTnv*v.*abs(n); % [N] thrust
Q_fit = cQn*n.*abs(n)-cQnv*v.*abs(n); % [Nm] torque

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