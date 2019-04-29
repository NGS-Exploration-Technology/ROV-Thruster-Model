clear all;
% Initialization of physical model constants
thruster_config;
rho = 1027;
Diam = Thruster_Config.D;
L = Thruster_Config.L;
kn = Thruster_Config.kn1;
kq = Thruster_Config.kq*2*pi/60;
kv1 = 7578*2*pi/60;
kv2 = 6770*2*pi/60;
dv1 = -.8475;
dv2 = .9254;
udotmax = 1/.01;
cT1 = Thruster_Config.cT1/((2*pi/60)^2);
cQ1 = Thruster_Config.cQ1/((2*pi/60)^2);

% Initialization of controller parameters
Kn = 10;
cntrl_dt = .05; % Sampling frequency
rate = robotics.Rate(1/cntrl_dt);
t = 0:cntrl_dt:4;
Td = 10*ones(size(t)); %*sin(0.5*t.');
nd = zeros(length(t),1);
nhat = zeros(length(t),1);
u = 0;
uprev = 0;
y_volts = zeros(length(t),1);

%%%%% NEW COVARIANCES FOR KALMAN FILTER %%%%%
Q = 1;
R = 2;
P = R;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize Arduino Communication
% UPDATE WITH CURRENT ARDUINO COM PORT
%ard = serial('COM4', 'BaudRate', 9600, 'DataBits', 8);
TACH = 'A0';
FLOW = 'A1';
OUT = 'D5';
ard = arduino('COM4','uno');
configurePin(ard,OUT,'PWM');
configurePin(ard,TACH,'AnalogInput');
configurePin(ard,FLOW,'AnalogInput');
configurePin(ard,'D6','DigitalOutput');
configurePin(ard,'D7','DigitalOutput');


%Set initial voltage to 0V
writePWMVoltage(ard,OUT,2.5);

% Enable output from Arduino
writeDigitalPin(ard,'D6',0);
writeDigitalPin(ard,'D7',1);

fprintf(1, 'Entering Control Loop!\n');
keyboard
reset(rate);

tic
for i = 1:length(Td)
    % read tachometer and flow speed
    y_volts(i) = readVoltage(ard,TACH);
    y = (477.43*(2*y_volts(i)-5)+5.3255)*2*pi/60;
    % flow = readVoltage(ard,FLOW);

    % Estimate update
    K = P/(P+R);
    nhat(i) = nhat(i)+K*(sign(nhat(i))*abs(y)-nhat(i));
    P = (1-K)*P;

    % Control update
    nd(i) = sqrt(abs(Td(i)/(rho*Diam^4*cT1)))*sign(Td(i));
    ndot_d = (nd(i)-nhat(i))/cntrl_dt;
    a = ndot_d+kn*nd(i)-kq*rho*Diam^5*cQ1*nd(i)*abs(nd(i));
    b = Kn*(nhat(i)-nd(i));
    if a < b
        u = ((a-b)/kv1)+dv1;
    elseif a > b
        u = ((a-b)/kv2)+dv2;
    end
    if (abs(u-uprev)/cntrl_dt)>udotmax
        u = uprev+sign(u-uprev)*udotmax*cntrl_dt;  
    end
    if abs(u) > 3
        u = 3*sign(u);
    end
    if u <= dv1
        gamma = kv1 * (u - dv1);
    elseif u >= dv2
        gamma = kv2 * (u - dv2);
    else
        gamma = 0;
    end
    uprev = u;

    % u -> volts -> state_value (0-255)
    %    -> shift u into 0-10V range
    v_out = (-u + 5.0) / 2.0;
    writePWMVoltage(ard,OUT,v_out);

    % Observer prediction
    nhatdot = -kn*nhat(i)+kq*rho*Diam^5*cQ1*nhat(i)*abs(nhat(i))+gamma;
    nhat(i+1) = nhat(i)+nhatdot*cntrl_dt;
    A = -kn+2*kq*rho*Diam^5*cQ1*abs(nhat(i));
    Pdot = 2*A*P+Q;
    P = P+Pdot*cntrl_dt;

    % Wait for next loop
    waitfor(rate);

end
toc
writePWMVoltage(ard,OUT,2.5);
statistics(rate)

writeDigitalPin(ard,'D6',0);
writeDigitalPin(ard,'D7',0);
