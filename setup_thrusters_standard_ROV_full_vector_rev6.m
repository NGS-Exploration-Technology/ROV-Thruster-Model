function [Thruster_Config] = setup_thrusters_standard_ROV_full_vector_rev6()
%SETUP_THRUSTERS Generate input data structures for setup_sim_ROV
%[Thruster_Config] = setup_thrusters_standard_ROV_full_vector()
%
%inputs:
%   none
%
%outputs:
%   Thruster_Params = data structure defining thruster configuration for setup_sim_ROV
%
%See also SETUP_SIM_ROV SIM_ROV

Throttle_to_n_Table = csvread('Throttle_to_n_Table.csv');
Table_Throttles = Throttle_to_n_Table(:,1);
Table_n = Throttle_to_n_Table(:,2);

%Initialize Variables
Thruster_Count = 6;
index = 0;
while(index<Thruster_Count)
    index = index + 1;
    Thruster_Config(index).position = zeros(3,1); %[m] Location vector of thruster
    Thruster_Config(index).yaw_angle = 0; %[deg] Thruster yaw angle from x
    Thruster_Config(index).pitch_angle = 0; %[deg] Thruster pitch angle from xy plane
    Thruster_Config(index).max_thrust = 0; %[N] maximum thruster output
    Thruster_Config(index).d_thrust_dt = 0; %[N/s] maximum allowable change in thrust
    Thruster_Config(index).kv = 0; %RPM Throttle rate parameter
    Thruster_Config(index).cT1 = 0; %Piecewise Alpha1 reverse coefficient
    Thruster_Config(index).cT2 = 0; %Piecewise Alpha1 forward coefficient
    Thruster_Config(index).dT1 = 0; %Piecewise Alpha1 reverse deadband size
    Thruster_Config(index).dT2 = 0; %Piecewise Alpha1 forward deadband
    Thruster_Config(index).cQ1 = 0; %Piecewise Beta1 reverse coeficient
    Thruster_Config(index).cQ2 = 0; %Piecewise Beta1 forward coeficient
    Thruster_Config(index).dQ1 = 0; %Piecewise Beta1 reverse deadband
    Thruster_Config(index).dQ2 = 0; %Piecewise Beta1 forward deadband
    Thruster_Config(index).g = 0; %[N/Throttle] throttle gain parameter
    Thruster_Config(index).D = 0; %[m] diameter of the prop
    Thruster_Config(index).alpha2 = 0; %Thrust coefficient 2 (1 is calculated)
    Thruster_Config(index).beta2 = 0; %Torque coefficient 2 (1 is calculated)
    Thruster_Config(index).RH_prop = 0; %1 = RH prop, 0 = LH prop
    Thruster_Config(index).Table_Throttles = Table_Throttles; %throttle to n Table_Throttles (applied to each thruster)
    Thruster_Config(index).Table_n = Table_n; %throttle to n Table_n (applied to each thruster)
end

%Define Thruster Configuration (index for Thruster_Params is Thruster ID)
    %Thruster 1
    Thruster_Config(1).position(1)= 0.24892; %[m] x Location of Thruster on frame WRT center of mass
    Thruster_Config(1).position(2) = -0.219202; %[m] y Location of Thruster on frame
    Thruster_Config(1).position(3) = -0.019812; %[m] z Location of Thruster on frame
    Thruster_Config(1).yaw_angle = 30; %[deg] Angle of Thruster from X
    Thruster_Config(1).pitch_angle = -30; %[deg] Angle of Thruster from XY plane
    Thruster_Config(1).max_thrust = 57.83; %[N] maximum thruster output for thruster 1
    Thruster_Config(1).d_thrust_dt = 1156; %[N/s] maximum allowable change in thrust for thruster 1
    Thruster_Config(1).g = 32.5; %[N/Throttle] throttle gain parameter for thruster 1
    Thruster_Config(1).D = 0.1151; %[m] diameter of the prop
    Thruster_Config(1).kv = 15.1307; %RPM Throttle rate parameter
    Thruster_Config(1).cT1 = 5.8078E-5; %Piecewise Alpha1 reverse coefficient
    Thruster_Config(1).cT2 = 5.0183E-5; %Piecewise Alpha1 forward coefficient
    Thruster_Config(1).dT1 = -3.4063; %Piecewise Alpha1 reverse deadband size
    Thruster_Config(1).dT2 = 1.4330; %Piecewise Alpha1 forward deadband
    Thruster_Config(1).cQ1 = -8.3319E-6; %Piecewise Beta1 reverse coeficient
    Thruster_Config(1).cQ2 = -1.7209E-5; %Piecewise Beta1 forward coeficient
    Thruster_Config(1).dQ1 = -2.7242; %Piecewise Beta1 reverse deadband
    Thruster_Config(1).dQ2 = 0.1696; %Piecewise Beta1 forward deadband
    Thruster_Config(1).alpha2 = 0.0011; %Thrust coefficient 2 (1 is calculated)
    Thruster_Config(1).beta2 = 0.0030; %Torque coefficient 2 (1 is calculated)
    Thruster_Config(1).RH_prop = 0; %1 = RH prop, 0 = LH prop
    %Thruster 2
    Thruster_Config(2).position(1) = 0.24892; %[m] x Location of Thruster on frame
    Thruster_Config(2).position(2) = 0.219202; %[m] y Location of Thruster on frame
    Thruster_Config(2).position(3) = -0.019812; %[m] z Location of Thruster on frame
    Thruster_Config(2).yaw_angle = -30; %[deg] Angle of Thruster from X
    Thruster_Config(2).pitch_angle = -30; %[deg] Angle of Thruster from XY plane
    Thruster_Config(2).max_thrust = 57.83; %[N] maximum thruster output for thruster 1
    Thruster_Config(2).d_thrust_dt = 1156; %[N/s] maximum allowable change in thrust for thruster 1
    Thruster_Config(2).g = 32.5; %[N/Throttle] throttle gain parameter for thruster 1
    Thruster_Config(2).D = 0.1151; %[m] diameter of the prop
    Thruster_Config(2).kv = 15.1307; %RPM Throttle rate parameter
    Thruster_Config(2).cT1 = 5.8078E-5; %Piecewise Alpha1 reverse coefficient
    Thruster_Config(2).cT2 = 5.0183E-5; %Piecewise Alpha1 forward coefficient
    Thruster_Config(2).dT1 = -3.4063; %Piecewise Alpha1 reverse deadband size
    Thruster_Config(2).dT2 = 1.4330; %Piecewise Alpha1 forward deadband
    Thruster_Config(2).cQ1 = -8.3319E-6; %Piecewise Beta1 reverse coeficient
    Thruster_Config(2).cQ2 = -1.7209E-5; %Piecewise Beta1 forward coeficient
    Thruster_Config(2).dQ1 = -2.7242; %Piecewise Beta1 reverse deadband
    Thruster_Config(2).dQ2 = 0.1696; %Piecewise Beta1 forward deadband
    Thruster_Config(2).alpha2 = 0.0011; %Thrust coefficient 2 (1 is calculated)
    Thruster_Config(2).beta2 = 0.0030; %Torque coefficient 2 (1 is calculated)
    Thruster_Config(2).RH_prop = 1; %1 = RH prop, 0 = LH prop
    %Thruster 3
    Thruster_Config(3).position(1) = -0.20828; %[m] x Location of Thruster on frame
    Thruster_Config(3).position(2) = 0.2413; %[m] y Location of Thruster on frame
    Thruster_Config(3).position(3) = -0.019812; %[m] z Location of Thruster on frame
    Thruster_Config(3).yaw_angle = 30; %[deg] Angle of Thruster from X
    Thruster_Config(3).pitch_angle = -30; %[deg] Angle of Thruster from XY plane
    Thruster_Config(3).max_thrust = 57.83; %[N] maximum thruster output for thruster 1
    Thruster_Config(3).d_thrust_dt = 1156; %[N/s] maximum allowable change in thrust for thruster 1
    Thruster_Config(3).g = 32.5; %[N/Throttle] throttle gain parameter for thruster 1
    Thruster_Config(3).D = 0.1151; %[m] diameter of the prop
    Thruster_Config(3).kv = 15.1307; %RPM Throttle rate parameter
    Thruster_Config(3).cT1 = 5.8078E-5; %Piecewise Alpha1 reverse coefficient
    Thruster_Config(3).cT2 = 5.0183E-5; %Piecewise Alpha1 forward coefficient
    Thruster_Config(3).dT1 = -3.4063; %Piecewise Alpha1 reverse deadband size
    Thruster_Config(3).dT2 = 1.4330; %Piecewise Alpha1 forward deadband
    Thruster_Config(3).cQ1 = -8.3319E-6; %Piecewise Beta1 reverse coeficient
    Thruster_Config(3).cQ2 = -1.7209E-5; %Piecewise Beta1 forward coeficient
    Thruster_Config(3).dQ1 = -2.7242; %Piecewise Beta1 reverse deadband
    Thruster_Config(3).dQ2 = 0.1696; %Piecewise Beta1 forward deadband
    Thruster_Config(3).alpha2 = 0.0011; %Thrust coefficient 2 (1 is calculated)
    Thruster_Config(3).beta2 = 0.0030; %Torque coefficient 2 (1 is calculated)
    Thruster_Config(3).RH_prop = 0; %1 = RH prop, 0 = LH prop
    %Thruster 4
    Thruster_Config(4).position(1) = -0.20828; %[m] x Location of Thruster on frame
    Thruster_Config(4).position(2) = -0.2413; %[m] y Location of Thruster on frame
    Thruster_Config(4).position(3) = -0.019812; %[m] z Location of Thruster on frame
    Thruster_Config(4).yaw_angle = -30; %[deg] Angle of Thruster from X
    Thruster_Config(4).pitch_angle = -30; %[deg] Angle of Thruster from XY plane
    Thruster_Config(4).max_thrust = 57.83; %[N] maximum thruster output for thruster 1
    Thruster_Config(4).d_thrust_dt = 1156; %[N/s] maximum allowable change in thrust for thruster 1
    Thruster_Config(4).g = 32.5; %[N/Throttle] throttle gain parameter for thruster 1
    Thruster_Config(4).D = 0.1151; %[m] diameter of the prop
    Thruster_Config(4).kv = 15.1307; %RPM Throttle rate parameter
    Thruster_Config(4).cT1 = 5.8078E-5; %Piecewise Alpha1 reverse coefficient
    Thruster_Config(4).cT2 = 5.0183E-5; %Piecewise Alpha1 forward coefficient
    Thruster_Config(4).dT1 = -3.4063; %Piecewise Alpha1 reverse deadband size
    Thruster_Config(4).dT2 = 1.4330; %Piecewise Alpha1 forward deadband
    Thruster_Config(4).cQ1 = -8.3319E-6; %Piecewise Beta1 reverse coeficient
    Thruster_Config(4).cQ2 = -1.7209E-5; %Piecewise Beta1 forward coeficient
    Thruster_Config(4).dQ1 = -2.7242; %Piecewise Beta1 reverse deadband
    Thruster_Config(4).dQ2 = 0.1696; %Piecewise Beta1 forward deadband
    Thruster_Config(4).alpha2 = 0.0011; %Thrust coefficient 2 (1 is calculated)
    Thruster_Config(4).beta2 = 0.0030; %Torque coefficient 2 (1 is calculated)
    Thruster_Config(4).RH_prop = 0; %1 = RH prop, 0 = LH prop
    %Thruster 5
    Thruster_Config(5).position(1) = 0; %[m] x Location of Thruster on frame
    Thruster_Config(5).position(2) = -0.1651; %[m] y Location of Thruster on frame
    Thruster_Config(5).position(3) = -0.05715; %[m] z Location of Thruster on frame
    Thruster_Config(5).yaw_angle = 90; %[deg] Angle of Thruster from X
    Thruster_Config(5).pitch_angle = -85; %[deg] Angle of Thruster from XY plane
    Thruster_Config(5).max_thrust = 57.83; %[N] maximum thruster output for thruster 1
    Thruster_Config(5).d_thrust_dt = 1156; %[N/s] maximum allowable change in thrust for thruster 1
    Thruster_Config(5).g = 32.5; %[N/Throttle] throttle gain parameter for thruster 1
    Thruster_Config(5).D = 0.1151; %[m] diameter of the prop
    Thruster_Config(5).kv = 15.1307; %RPM Throttle rate parameter
    Thruster_Config(5).cT1 = 5.8078E-5; %Piecewise Alpha1 reverse coefficient
    Thruster_Config(5).cT2 = 5.0183E-5; %Piecewise Alpha1 forward coefficient
    Thruster_Config(5).dT1 = -3.4063; %Piecewise Alpha1 reverse deadband size
    Thruster_Config(5).dT2 = 1.4330; %Piecewise Alpha1 forward deadband
    Thruster_Config(5).cQ1 = -8.3319E-6; %Piecewise Beta1 reverse coeficient
    Thruster_Config(5).cQ2 = -1.7209E-5; %Piecewise Beta1 forward coeficient
    Thruster_Config(5).dQ1 = -2.7242; %Piecewise Beta1 reverse deadband
    Thruster_Config(5).dQ2 = 0.1696; %Piecewise Beta1 forward deadband
    Thruster_Config(5).alpha2 = 0.0011; %Thrust coefficient 2 (1 is calculated)
    Thruster_Config(5).beta2 = 0.0030; %Torque coefficient 2 (1 is calculated)
    Thruster_Config(5).RH_prop = 0; %1 = RH prop, 0 = LH prop
    %Thruster 6
    Thruster_Config(6).position(1) = 0; %[m] x Location of Thruster on frame
    Thruster_Config(6).position(2) = 0.1651; %[m] y Location of Thruster on frame
    Thruster_Config(6).position(3) = -0.05715; %[m] z Location of Thruster on frame
    Thruster_Config(6).yaw_angle = -90; %[deg] yaw angle of Thruster from X
    Thruster_Config(6).pitch_angle = -85; %[deg] pitch ngle of Thruster from XY plane
    Thruster_Config(6).max_thrust = 57.83; %[N] maximum thruster output for thruster 1
    Thruster_Config(6).d_thrust_dt = 1156; %[N/s] maximum allowable change in thrust for thruster 1
    Thruster_Config(6).g = 32.5; %[N/Throttle] throttle gain parameter for thruster 1
    Thruster_Config(6).D = 0.1151; %[m] diameter of the prop
    Thruster_Config(6).kv = 15.1307; %RPM Throttle rate parameter
    Thruster_Config(6).cT1 = 5.8078E-5; %Piecewise Alpha1 reverse coefficient
    Thruster_Config(6).cT2 = 5.0183E-5; %Piecewise Alpha1 forward coefficient
    Thruster_Config(6).dT1 = -3.4063; %Piecewise Alpha1 reverse deadband size
    Thruster_Config(6).dT2 = 1.4330; %Piecewise Alpha1 forward deadband
    Thruster_Config(6).cQ1 = -8.3319E-6; %Piecewise Beta1 reverse coeficient
    Thruster_Config(6).cQ2 = -1.7209E-5; %Piecewise Beta1 forward coeficient
    Thruster_Config(6).dQ1 = -2.7242; %Piecewise Beta1 reverse deadband
    Thruster_Config(6).dQ2 = 0.1696; %Piecewise Beta1 forward deadband
    Thruster_Config(6).alpha2 = 0.0011; %Thrust coefficient 2 (1 is calculated)
    Thruster_Config(6).beta2 = 0.0030; %Torque coefficient 2 (1 is calculated)
    Thruster_Config(6).RH_prop = 1; %1 = RH prop, 0 = LH prop
    
    %Calculate unit direction vectors for each thruster
    Thruster_Count = length(Thruster_Config);
    u_hat = zeros(3,1);
    index = 0;
    while(index<Thruster_Count)
        index = index + 1;
        u_hat(1) = cos(Thruster_Config(index).yaw_angle*pi/180)*cos(Thruster_Config(index).pitch_angle*pi/180); %Unit vector x component (calculated from angles)
        u_hat(2) = sin(Thruster_Config(index).yaw_angle*pi/180)*cos(Thruster_Config(index).pitch_angle*pi/180); %Unit vector y component (calculated from angles)
        u_hat(3) = sin(Thruster_Config(index).pitch_angle*pi/180); %Unit vector z component (calculated from pitch angles)
        Thruster_Config(index).direction = u_hat;
    end
    