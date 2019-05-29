clear

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
Kp = 10;
Ki = 10;
dt = .03; % Sampling frequency
rate = robotics.Rate(1/dt);
t = 0:dt:10;
% Td = [10*ones(length(t)/2,1);-10*ones(length(t)/2,1)];
% Td = [(linspace(0,10,50)).';(linspace(10,-10,100)).';(linspace(-10,0,50)).';zeros(51,1)];
Td = 10*sin(1.88*t.');
nd = zeros(length(t),1);
nhat = zeros(length(t),1);
u = zeros(length(t),1);
uprev = 0;
udotmax = 1/.01;
% load('y_volts_step_updown.mat');
load('y_volts_sinusoid.mat');
y = 49.9723*(2*y_volts-5);
c = 0;

% Covariances for KF:
Q = 1;
R = 2;
P = R;

% Initialize Arduino Communication
% UPDATE WITH CURRENT ARDUINO COM PORT
% ard = serial('COM4', 'BaudRate', 9600, 'DataBits', 8);
% TACH = 'A0';
% FLOW = 'A1';
% OUT = 'D5';
% ard = arduino('COM4','uno');
% configurePin(ard,OUT,'PWM');
% configurePin(ard,TACH,'AnalogInput');
% configurePin(ard,FLOW,'AnalogInput');
% configurePin(ard,'D6','DigitalOutput');
% configurePin(ard,'D7','DigitalOutput');
% %Set initial voltage to 0V
% writePWMVoltage(ard,OUT,2.5);
% % Enable output from Arduino
% writeDigitalPin(ard,'D6',0);
% writeDigitalPin(ard,'D7',1);
% 
% fprintf(1, 'Entering Control Loop!\n');
% s = serial('COM4', 'BaudRate', 9600);
% fopen(s);
% 
% keyboard
% reset(rate);
% tic
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
    
    % u -> volts -> state_value (0-255)
    %    -> shift u into 0-10V range
    % v_out = (-u(i) + 5.0);
    % v_out = v_out * 255 / 10;
    % writePWMVoltage(ard,OUT,v_out);
    % fprintf(s, '%d\n', round(v_out));
    
    % read tachometer and flow speed
    % y_volts(i) = readVoltage(ard,TACH);
    % flow = readVoltage(ard,FLOW);
    % reading = fscanf(s, '%d,%d');
    % y_raw = reading(1);
    % flow_raw = reading(2);
    % y_volts = y_raw * 5 / 1024;
    % flow_volts = flow_raw * 5 / 1024;
    % y(i) = 49.9723*(2*y_volts-5);
    
    % Observer prediction
    K = P/R;
    nhatdot = -kn*nhat(i)-kq*cQn*nhat(i)*abs(nhat(i))+gamma+K*(sign(nhat(i))*abs(y(i))-nhat(i));
    nhat(i+1) = nhat(i)+nhatdot*dt;
    A = -kn-2*kq*cQn*abs(nhat(i));
    Pdot = 2*A*P-K*P+Q;
    P = P+Pdot*dt;

    % Wait for next loop
    % waitfor(rate);
end
% toc
% statistics(rate)
% writePWMVoltage(ard,OUT,2.5);
% reading = fscanf(s, '%d,%d');
% fprintf(s, '-99\n');
% fclose(s);

% writeDigitalPin(ard,'D6',0);
% writeDigitalPin(ard,'D7',0);
figure
plot(t,y,t,abs(nhat(1:(end-1))),t,abs(nd),'--k')
axis([0 10 0 110])
xlabel('Time [s]')
ylabel('Propeller Speed [rad/s]')
legend({'Tachometer Reading, $|n|$','Estimate Magnitude, $|\hat{n}|$','Setpoint Magnitude, $|n_{d}|$'},'Interpreter','Latex')
grid
figure
plot(t,y,t,nhat(1:(end-1)),t,nd,'--k')
axis([0 10 -110 110])
xlabel('Time [s]')
ylabel('Propeller Velocity [rad/s]')
legend({'Tachometer Reading, n','Estimate, $\hat{n}$','Setpoint, $n_{d}$'},'Interpreter','Latex')
grid
% figure
% plot(t,u)
% grid