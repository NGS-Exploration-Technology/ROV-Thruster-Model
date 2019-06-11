clear; close all;

% Initialization of physical model constants
thruster_config;
cTn = Thruster_Config.cTn;
cQn = Thruster_Config.cQn;

% Initialization of controller parameters
dt = .04; % Sampling frequency
rate = robotics.Rate(1/dt);
t = 0:dt:10;
% Td = -[(linspace(0,10,50)).';10*ones(((length(t)-1)/2)-100,1);(linspace(10,-10,100)).';-10*ones(((length(t)+1)/2)-100,1);(linspace(-10,0,50)).'];
% Td = [-10*ones((length(t)-1)/2,1);10*ones((length(t)+1)/2,1)];
Td = 5*cos(1.88*t.')-5;
nd = zeros(length(t),1);
u = zeros(length(t),1);
uprev = 0;
udotmax = 1/.01;
y = zeros(length(t),1);

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
    if nd(i) < 0
        u(i) = (nd(i)-38.756)/50.014;
    elseif nd(i) > 0
        u(i) = (nd(i)+40.4)/50.031;
    end
    if (abs(u(i)-uprev)/dt)>udotmax
        u(i) = uprev+sign(u(i)-uprev)*udotmax*dt;  
    end
    if abs(u(i)) > 3
        u(i) = 3*sign(u(i));
    end
    uprev = u(i);
    
    % u -> volts -> state_value (0-255) -> shift u into 0-10V range
    v_out = (-u(i) + 5.0);
    v_out = v_out * 255 / 10;
    fprintf(s, '%d\n', round(v_out));
    
    % read tachometer and flow speed
    reading = fscanf(s, '%d,%d');
    y_raw = reading(1);
    y(i) = 49.9723*y_raw*5/1024;
    
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
subplot(2,1,1)
plot(t,y,t,abs(nd),'--k')
ylabel('Propeller Speed [rad/s]')
legend({'y','$|n_{d}|$'},'Interpreter','Latex')
grid
subplot(2,1,2)
plot(t,u)
ylabel('Input [V]')
xlabel('Time [s]')
grid