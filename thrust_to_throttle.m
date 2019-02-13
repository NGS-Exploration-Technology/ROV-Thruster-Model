function Throttle=thrust_to_throttle(Thrust, Table_Thrusts, Table_Throttles)
%function Throttle=thrust_to_throttle(Thrust, Table_Thrusts, Table_Throttles)
%THRUST_TO_THROTTLE calculate the Throttle based on a desired Thrust
%
%Inputs:
%   Thrust = [N] Thruster Output
%   Thrust_Table = [N] Vector of Thrusts associated with Throttle Table
%   Throttle_Table = Vector of Throttles associated with Thrust Table
%
%Outputs
%   Throttle = Throttle Value from -1 to 1
%
%National Geographic Society
%February 13, 2019
%Eric Berkenpas

min_Thrusts = Table_Thrusts(1);
max_Thrusts = Table_Thrusts(end);

%Pad Table_Throttles with Last Value
Table_Throttles(end+1) = Table_Throttles(end);

%Calculate fractional table index
Table_Index = ((Thrust-min_Thrusts)/(max_Thrusts-min_Thrusts))*127+1;

%Saturation Filter
if(Table_Index>128)
    Table_Index = 128;
end
if(Table_Index<1)
    Table_Index = 1;
end

%Calculate Real
Real = floor(Table_Index);

%Calculate Fraction
Fraction = Table_Index-Real;

%Get First AD Cts
Table_Point1 = Table_Throttles(Real);

%Get Second AD Cts
Table_Point2 = Table_Throttles(Real+1);

%Interpolate Between Fraction
Interpolated_Throttle = (Table_Point2-Table_Point1)*Fraction;

%Add Interpolated Counts to First AD Cts
Throttle = Table_Point1+Interpolated_Throttle;

%Test Point (debugging)
%Fbc = -sign(Velocity).*rho*A*Cd.*Velocity.^2;
%Position_actual = Fbc.*steps_per_newton;
