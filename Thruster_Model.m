function [n, T, Q, Va] = Thruster_Model(Va, n0, T0, Q0, Throttle , rho, Thruster_Config, dt)
%THRUSTER_MODEL time varying dynamic model of ROV thruster
%function [n, T, Q] = Thruster_Model(Va, n0, Throttle , rho, Thruster_Params, dt)
% T(t)=Ta+(To-Ta)*exp(-kt) 
%
%inputs:
%   Va = [m/s] Flow velocity at prop
%   n0 = [rps] current revolutions per second
%   T0 = [N] current thrust
%   Q0 = [Nm] current torque
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
% g = Thruster_Config.g; %[rps/throttle]
% k = Thruster_Config.k; %[rate parameter]
D = Thruster_Config.D; %[m] propellor diameter
kt = Thruster_Config.kt;
kn1 = Thruster_Config.kn1; %[rate parameter]
kn2 = Thruster_Config.kn2; %[rate parameter]
kq = Thruster_Config.kq; %[rate parameter]
kv = Thruster_Config.kv; %[rate parameter]
ku1 = Thruster_Config.ku1; %[rate parameter]
ku2 = Thruster_Config.ku2; %[rate parameter]
cT1 = Thruster_Config.cT1;
cT2 = Thruster_Config.cT2;
dT1 = Thruster_Config.dT1;
dT2 = Thruster_Config.dT2;
cQ1 = Thruster_Config.cQ1;
cQ2 = Thruster_Config.cQ2;
dQ1 = Thruster_Config.dQ1;
dQ2 = Thruster_Config.dQ2;
alpha2 = Thruster_Config.alpha2;
beta2 = Thruster_Config.beta2;

%Update thruster rate
% n_command = g*Throttle;
% d_n = -k*(n0-n_command);
if Thruster_Config.RH_prop
    d_n = -kn1*n0-kn2*n0*abs(n0)+kq*Q0+kv*Throttle; % RH prop
else
    d_n = -kn1*n0-kn2*n0*abs(n0)-kq*Q0+kv*Throttle; % LH prop
end
n = n0+d_n*dt;
d_u = -ku1*Va-ku2*Va*abs(Va)+kt*T0;
Va = Va+d_u*dt;

%Calculate Advance Coefficient
if (n==0) 
   n = 1E-100;
end
J0 = Va/(n*D);
%J0 = 0; %Assume the flow velocity is zero

% Calculate deadband piecewise quadratic coefficients:
if n*abs(n)<=dT1
    alpha1 = cT1*((n*abs(n))-dT1);
elseif n*abs(n)>=dT2
    alpha1 = cT2*((n*abs(n))-dT2);
else
    alpha1 = 0;
end
if n*abs(n)<=dQ1
    beta1 = cQ1*((n*abs(n))-dQ1);
elseif n*abs(n)>=dQ2
    beta1 = cQ2*((n*abs(n))-dQ2);
else
    beta1 = 0;
end

%Calculate Thrust
T = rho*D^4*(alpha1-alpha2*J0*abs(n)*n); %[N] thrust

%Calculate Torque
if Thruster_Config.RH_prop
    Q = rho*D^5*(beta1+beta2*J0*abs(n)*n); %[Nm] torque
else
    Q = rho*D^5*(beta1-beta2*J0*abs(n)*n); %[Nm] torque
end
end
