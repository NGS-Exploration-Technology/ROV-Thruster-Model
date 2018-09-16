clear; close all;
u_d = -5:.5:5;
rpm_d = [-2305.4357 -2071.495 -1737.294 -1431.7388 -1150.0551 -892.2429 -634.4307 -371.8442 -118.8063 0 0 0 105.5858 358.6237 606.8873 888.571 1146.3832 1432.8412 1747.945 2129.889 2425.8956];

u_c = -5:.001:5;
rpm_c = [554.02*u_c(u_c<=(-468.7/553.02))+468.7,0*u_c(u_c>(-468.7/553.02) & u_c<(537.41/580.71)),580.71*u_c(u_c>=(537.41/580.71))-537.41];

font = 14;
figure
plot(u_d,rpm_d,'o',u_c,rpm_c,'--k')
xlabel('Throttle [V]','FontSize',font,'FontName','Times New Roman')
xticks(-5:5)
ylabel('Propeller Speed [rpm]','FontSize',font,'FontName','Times New Roman')
legend({'Experimental','Model'},'FontSize',font,'FontName','Times New Roman','Location','Northwest')
set(gca,'FontSize',font,'FontName','Times New Roman');
grid