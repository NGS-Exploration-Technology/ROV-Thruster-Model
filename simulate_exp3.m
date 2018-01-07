%Approximate (simulate) an exponential function by integrating the derivative
%
%function to approximate T(t)=Ta+(To-Ta)*exp^{-kt}  
%
%dT/dt=-k(T-Ta). 
%T(0)=To  
%
%http://www.ugrad.math.ubc.ca/coursedoc/math100/notes/diffeqs/cool.html

%Simulation
dt = 0.1;
t_n = [0:dt:5];
To = 0;
Ta = 1;
k = 1;
index = 1;
T_n(index) = To;
while(t_n(index)<max(t_n))
    dT(index) = -k*(T_n(index)-Ta);
    T_n(index+1) = T_n(index)+dT(index)*dt;
    index = index + 1;
end

%Calculate true Exponential
t = linspace(0,5,1000);
T=Ta+(To-Ta)*exp(-k*t);

%Plot results
plot(t, T, t_n, T_n, '.k');