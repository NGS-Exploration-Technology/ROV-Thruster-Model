function [T, tau] = sim_thruster(To, Ta , k, g, I, v, dt)
%SIM_THRUSTER to fit a thrust throttle curve to a model
% T(t)=Ta+(To-Ta)*exp(-kt) 
%
%inputs:
%   To = [N] Current thrust
%   Ta = [N] Thrust Commanded (Throttle) Vector [N]
%   k = rate parameter
%   g = [N/Throttle] throttle gain parameter
%   I = Moment of inertia for the rotating machinery
%   Cd = Coefficient of rotational drag as a function of thrust
%   v = current flow velocity over the thruster
%   dt = simulation time step
%
%outputs:
%   T = [N] Thrust
%   tau = [Nm] Torque
%
%function to approximate T(t)=Ta+(To-Ta)*exp^{-kt}  
%
%dT/dt=-k(T-Ta). 
%T(0)=To  
%
%http://www.ugrad.math.ubc.ca/coursedoc/math100/notes/diffeqs/cool.html

dT = -k*(To-g*Ta);
T = To+dT*dt;
tau_inertial = I*(dT/dt);
tau_viscous = Cd*T;
tau = tau_inertial+tau_viscous;