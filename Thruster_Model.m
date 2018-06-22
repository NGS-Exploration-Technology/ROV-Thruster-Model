function [n, T, Q] = Thruster_Model(Va, n0, Throttle , rho, Thruster_Config, dt)
%THRUSTER_MODEL time varying dynamic model of ROV thruster
%function [n, T, Q] = Thruster_Model(Va, n0, Throttle , rho, Thruster_Params, dt)
% T(t)=Ta+(To-Ta)*exp(-kt) 
%
%inputs:
%   Va = [m/s] Flow velocity at prop
%   n0 = [rps] current revolutions per second
%   Throttle = current throttle (-1 to 1)
%   rho = density of seawater
%   Thruster_Params = parameters of the thruster model
%   dt = simulation time step
%
%outputs:
%   n = [rps] thruster rate
%   T = [N] Output Thrust
%   Q = [Nm] Reaction Torque
%
%function to approximate T(t)=Ta+(To-Ta)*exp^{-kt}  
%
%dT/dt=-k(T-Ta). 
%T(0)=To  
%
%http://www.ugrad.math.ubc.ca/coursedoc/math100/notes/diffeqs/cool.html

%Populate thruster params
g = Thruster_Config.g; %[rps/throttle]
%k = Thruster_Config.k; %[rate parameter]
k1 = Thruster_Config.k1; %[rate parameter]
k2 = Thruster_Config.k2; %[rate parameter]
kt = Thruster_Config.kt; %[rate parameter]
D = Thruster_Config.D; %[m] propellor diameter
alpha1 = Thruster_Config.alpha1;
alpha2 = Thruster_Config.alpha2;
beta1 = Thruster_Config.beta1;
beta2 = Thruster_Config.beta2;

%Update thruster rate
n_command = g*Throttle;
% d_n = -k*(n0-n_command);
d_n = -k1*n0-k2*n0*abs(n0)-kt*n_command;
n = n0+d_n*dt;

%Calculate Advance Coefficient
%if (n==0) 
%    n = 1E-100;
%end
%J0 = Va/(n*D);
J0 = 0; %Assume the flow velocity is zero

%Calculate Thrust
T = rho*D^4*(alpha1+alpha2*J0)*abs(n)*n; %[N] thrust

%Calculate Torque
if Thruster_Config.RH_prop
    Q = rho*D^5*(beta1+beta2*J0)*abs(n)*n; %[Nm] RH prop torque
else
    Q = -rho*D^5*(beta1+beta2*J0)*abs(n)*n; %[Nm] LH prop torque
end
