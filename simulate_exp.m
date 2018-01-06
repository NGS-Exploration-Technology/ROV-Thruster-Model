%Approximate (simulate) an exponential function by integrating the derivative
%
%function to approximate y(t) = y0*exp(k*t)

%Simulation
dt = 0.1;
t_n = [0:dt:5];
k = -1;
y0 = 1;
index = 1;
y_n(index) = y0;
while(t_n(index)<max(t_n));
    dy(index) = k*y_n(index)*dt;
    y_n(index+1) = y_n(index)+dy(index);
    index = index + 1;
end

%Calculate true Exponential
t = linspace(0,5,1000);
y = y0*exp(k*t);

%Plot results
plot(t, y, t_n, y_n, '.k');