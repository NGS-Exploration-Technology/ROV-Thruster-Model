function [t_n,T_out] = sim_thruster(To, Ta , k, g, I, Cd, v, dt)
%SIM_THRUSTER to fit a thrust throttle curve to a model
% T(t)=Ta+(To-Ta)*exp(-kt) 
%
%inputs:
%   To = [N] Initial thrust
%   Ta = [N] Thrust Commanded (Throttle) Vector [N]
%   k = rate parameter
%   g = [N/Throttle] throttle gain parameter
%   I = Moment of inertia for the rotating machinery
%   Cd = Coefficient of rotational drag as a function of thrust
%   v = current flow velocity over the thruster
%   dt = simulation time step
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
    [T_out(index+1), tau(index+1)] = Thruster_Model(T_out(index), Ta(index), k, g, I, Cd, v, dt);
    index = index + 1;
end
