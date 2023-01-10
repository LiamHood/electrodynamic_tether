clear; close all; clc;
a0 = 6878;
e0 = .02;
i0_deg = 30;
i0 = deg2rad(i0_deg);
RAAN0_deg = 30;
RAAN0 = deg2rad(RAAN0_deg);
aop0_deg = 50;
aop0 = deg2rad(aop0_deg);
ta0_deg = 0;
ta0 = deg2rad(ta0_deg);

tend = 24*10;
mu = 396800;
forces = "gravityJ2J3";
dt = 100;
A = 1;
m = 1;

[ r0 , v0 ] = coes2state( [sqrt(mu*a0*(1-e0^2)), i0, e0, RAAN0, aop0, ta0], mu );
tspan = [0, tend*60*60];

[ t , rvec, vvec ] = Encke( dt , tspan , r0 , v0 , mu , forces , A , m );
t = t./(60*60);
for ii = 1:length(t)
    COES = state2coes_display( [rvec(:,ii), vvec(:,ii)], mu );
    a(ii) = COES(9);
    e(ii) = COES(3);
    i(ii) = COES(2);
    RAAN(ii) = COES(4);
    aop(ii) = COES(5);
    ta(ii) = COES(6);
    r(ii) = norm(rvec(:,ii));
    v(ii) = norm(vvec(:,ii));
end

figure
plot3(rvec(1,:), rvec(2,:), rvec(3,:))
xlabel("x")
ylabel("y")
zlabel("z")
axis equal

figure
subplot(2,1,1)
plot(t,r)
xlabel('time [hours]')
ylabel('Radius')

subplot(2,1,2)
plot(t,v)
xlabel('time [hours]')
ylabel('Velocity')

figure
subplot(3,2,1)
plot(t, a)
xlabel('time [hours]')
ylabel('Semi-Major Axis [km]')

subplot(3,2,2)
plot(t, e)
xlabel('time [hours]')
ylabel('Eccentricity')

subplot(3,2,3)
plot(t, i)
xlabel('time [hours]')
ylabel('Inclination [degree]')

subplot(3,2,4)
plot(t, RAAN)
xlabel('time [hours]')
ylabel('Right Ascension of Ascending Node [degree]')

subplot(3,2,5)
plot(t, aop)
xlabel('time [hours]')
ylabel('Argument of Perigee [degree]')

subplot(3,2,6)
plot(t, ta)
xlabel('time [hours]')
ylabel('True Anomaly [degree]')