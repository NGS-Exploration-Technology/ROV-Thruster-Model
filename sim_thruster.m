function [t_n, T, Q, n] = sim_thruster(n_command , rho, Thruster_Params, dt)
%SIM_THRUSTER to fit a thrust throttle curve to a model
%function [t_n,T, Q] = sim_thruster(Va, Throttle , rho, Thruster_Params, dt)
% T(t)=Ta+(To-Ta)*exp(-kt) 
%
%inputs:
%   Va = flow velocity at prop
%   Throttle = throttle setting (-1 to +1)
%   rho = density of fluid
%   Thruster_Params = data structure of params (g, k, D, alpha1, alpha2, beta1, beta2)
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

n_samples = length(n_command);
t_n = 0:dt:(n_samples*dt-dt);
n = zeros(size(t_n));
T = zeros(size(t_n));
Q = zeros(size(t_n));
index = 1;
while(t_n(index)<max(t_n))
    [n(index+1), T(index+1), Q(index+1)] = Thruster_Model(n(index), T(index), Q(index), n_command(index), rho, Thruster_Params, dt);
    index = index + 1;
end
