function [Thr, tau] = Thruster_Model(To, Ta , k, g, I, Cd, v, dt)
%THRUSTER_MODEL time varying dynamic model of ROV thruster
%function [T, tau] = Thruster_Model(To, Ta , k, g, I, Cd, v, dt)
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
%   Thr = [N] Output Thrust
%   tau = [Nm] Reaction Torque
%
%function to approximate T(t)=Ta+(To-Ta)*exp^{-kt}  
%
%dT/dt=-k(T-Ta). 
%T(0)=To  
%
%http://www.ugrad.math.ubc.ca/coursedoc/math100/notes/diffeqs/cool.html

dThr = -k*(To-g*Ta);
Thr = To+dThr*dt;
tau_inertial = I*(dThr/dt);
tau_viscous = Cd*Thr;
tau = tau_inertial+tau_viscous;