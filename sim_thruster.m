function [t_n,T_out] = sim_thruster(To, Ta , k, dt)
%SIM_THRUSTER to fit a thrust throttle curve to a model
% T(t)=Ta+(To-Ta)*exp(-kt) 
%
%inputs:
%   To = [N] Initial thrust
%   Ta = [N] Thrust Commanded (Throttle) Vector [N]
%   k = rate parameter
%   n = 
%
%outputs:
%
%function to approximate T(t)=Ta+(To-Ta)*exp^{-kt}  
%
%dT/dt=-k(T-Ta). 
%T(0)=To  
%
%http://www.ugrad.math.ubc.ca/coursedoc/math100/notes/diffeqs/cool.html

%set parameters

n = length(Ta);
t_n = [0:dt:n*dt];
T_out = zeros(size(t_n));
index = 1;
T_out(index) = To;
while(t_n(index)<max(t_n))
    dT(index) = -k*(T_out(index)-Ta(index));
    T_out(index+1) = T_out(index)+dT(index)*dt;
    index = index + 1;
end
