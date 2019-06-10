clear; close all;

% Initialization of physical model constants
thruster_config;
kn = Thruster_Config.kn;
kq = Thruster_Config.kq;
kv1 = Thruster_Config.kv1;
kv2 = Thruster_Config.kv2;
dv1 = Thruster_Config.dv1;
dv2 = Thruster_Config.dv2;
cTn = Thruster_Config.cTn;
cQn = Thruster_Config.cQn;

% Initialization of controller parameters
Kp = 5;
Ki = 10;
dt = .04; % Sampling frequency
rate = robotics.Rate(1/dt);
t = 0:dt:10;
% Td = -[(linspace(0,10,50)).';10*ones(((length(t)-1)/2)-100,1);(linspace(10,-10,100)).';-10*ones(((length(t)+1)/2)-100,1);(linspace(-10,0,50)).'];
Td = [-10*ones((length(t)-1)/2,1);10*ones((length(t)+1)/2,1)];
% Td = 5*cos(1.88*t.')-5;
nd = zeros(length(t),1);
nhat = zeros(length(t)+1,1);
u = zeros(length(t),1);
uprev = 0;
udotmax = 1/.01;
y = zeros(length(t),1);
c = 0;

% Covariances for KF:
QQ = 300;
RR = 3;
P = [RR;zeros(length(t),1)];

% Initialize Arduino Communication
fprintf(1, 'Entering Control Loop!\n');
s = serial('COM8', 'BaudRate', 9600);
fopen(s);
keyboard
reset(rate);
tic

% Enter loop
for i = 1:length(Td)
    
    % Control update
    nd(i) = sqrt(abs(Td(i)/cTn))*sign(Td(i));
    ndot_d = (nd(i)-nhat(i))/dt;
    cdot = nhat(i) - nd(i);
    c = c + cdot*dt;
    a = ndot_d+kn*nd(i)+kq*cQn*nd(i)*abs(nd(i));
    b = Kp*(nhat(i)-nd(i))+Ki*c;
    if a < b
        u(i) = ((a-b)/kv1)+dv1;
    elseif a > b
        u(i) = ((a-b)/kv2)+dv2;
    end
    if (abs(u(i)-uprev)/dt)>udotmax
        u(i) = uprev+sign(u(i)-uprev)*udotmax*dt;  
    end
    if abs(u(i)) > 3
        u(i) = 3*sign(u(i));
    end
    if u(i) <= dv1
        gamma = kv1*(u(i)-dv1);
    elseif u(i) >= dv2
        gamma = kv2*(u(i)-dv2);
    else
        gamma = 0;
    end
    uprev = u(i);
    
    % u -> volts -> state_value (0-255) -> shift u into 0-10V range
    v_out = (-u(i)+5)*255/10;
    fprintf(s, '%d\n', round(v_out));
    
    % read tachometer and flow speed
    reading = fscanf(s, '%d,%d');
    y(i) = 49.9723*reading(1)*5/1024;
    
    % Observer prediction
    K = P(i)/RR;
    nhatdot = -kn*nhat(i)-kq*cQn*nhat(i)*abs(nhat(i))+gamma+K*(sign(nhat(i))*abs(y(i))-nhat(i));
    nhat(i+1) = nhat(i)+nhatdot*dt;
    A = -kn-2*kq*cQn*abs(nhat(i));
    Pdot = 2*A*P(i)-K*P(i)+QQ;
    P(i+1) = P(i)+Pdot*dt;
    
    % Wait for next loop
    waitfor(rate);
end

toc
statistics(rate)
fprintf(s, '-99\n');
reading = fscanf(s, '%d,%d');
fclose(s);

% Plot results
figure
subplot(3,1,1)
plot(t,y,t,abs(nhat(1:(end-1))),t,abs(nd),'--k')
ylabel('Propeller Speed [rad/s]')
legend({'y','$|\hat{n}|$','$|n_{d}|$'},'Interpreter','Latex')
grid
subplot(3,1,2)
plot(t,P(1:(end-1)),'g')
ylabel('Variance')
grid
subplot(3,1,3)
plot(t,u)
ylabel('Input [V]')
xlabel('Time [s]')
grid