function n=throttle_to_n(Throttle, Table_Throttles, Table_n)
%THRUST_TO_THROTTLE calculate the Throttle based on a desired Thrust
%function n=throttle_to_n(Throttle, Table_Throttles, Table_n)
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

min_Throttles = Table_Throttles(1);
max_Throttles = Table_Throttles(end);

%Pad Table_Throttles with Last Value
Table_n(end+1) = Table_n(end);

%Calculate fractional table index
Table_Index = ((Throttle-min_Throttles)/(max_Throttles-min_Throttles))*127+1;

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
Table_Point1 = Table_n(Real);

%Get Second AD Cts
Table_Point2 = Table_n(Real+1);

%Interpolate Between Fraction
Interpolated_n = (Table_Point2-Table_Point1)*Fraction;

%Add Interpolated Counts to First AD Cts
n = Table_Point1+Interpolated_n;

%Test Point (debugging)
%Fbc = -sign(Velocity).*rho*A*Cd.*Velocity.^2;
%Position_actual = Fbc.*steps_per_newton;
