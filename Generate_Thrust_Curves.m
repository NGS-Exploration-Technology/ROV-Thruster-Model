function [cTn,cQn] = Generate_Thrust_Curves()

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
n = A(:,4);
T = A(:,2);
Q = A(:,3);

% fminsearch operation:
f = @(x)sqrt(mean((x(1)*n.*abs(n)-T).^2)+mean((x(2)*n.*abs(n)-Q).^2));

[x,fval,exitflag,opt_output] = fminsearch(f,[1 1], options)

cTn = x(1);
cQn = x(2);

% Calculate Thrust/Torque:
T_fit = cTn*n.*abs(n); % [N] thrust
Q_fit = cQn*n.*abs(n); % [Nm] torque

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