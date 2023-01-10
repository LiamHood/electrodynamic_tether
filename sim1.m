clear; close all; clc;

tspan = [0, 1e4];

a0 = 6878;
e0 = .02;
i0 = deg2rad(30);
RAAN0 = deg2rad(30);
aop0 = deg2rad(50);
ta0 = 0;
sc_state0 = [a0; e0; i0; RAAN0; aop0; ta0];

tether_state0 = [0.1; 0; 0; 0; 1];

L = 15;
m1 = 250;
m2 = 150;
mt = 20;
tether_param = [L; m1; m2; mt];

mu = 398600;
tol = 1e-10;
[ t , states] = BasicNaiveTether( tspan , sc_state0, tether_state0, ...
    tether_param, mu , tol );

t = t/(60*60);
a = states(:,1);
e = states(:,2);
i = rad2deg(states(:,3));
RAAN = rad2deg(states(:,4));
aop = rad2deg(states(:,5));
ta = rad2deg(states(:,6));

theta = rad2deg(states(:,7));
phi = rad2deg(states(:,8));
dtheta = rad2deg(states(:,9));
dphi = rad2deg(states(:,10));

figure
plot(t, a)
xlabel('time [hours]')
ylabel('Semi-Major Axis [km]')

figure
plot(t, e)
xlabel('time [hours]')
ylabel('Eccentricity')

figure
plot(t, i)
xlabel('time [hours]')
ylabel('Inclination [degree]')

figure
plot(t, RAAN)
xlabel('time [hours]')
ylabel('Right Ascension of Ascending Node [degree]')

figure
plot(t, aop)
xlabel('time [hours]')
ylabel('Argument of Perigee [degree]')

figure
plot(t, ta)
xlabel('time [hours]')
ylabel('True Anomaly [degree]')

figure
plot(t, theta)
xlabel('time [hours]')
ylabel('In Plane Libration [degree]')

figure
plot(t, phi)
xlabel('time [hours]')
ylabel('Out of Plane Libration [degree]')