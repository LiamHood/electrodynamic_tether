clear; close all; clc;

%% Defaults
a0 = 6978;
e0 = .02;
i0_deg = 30;
i0 = deg2rad(i0_deg);
RAAN0_deg = 30;
RAAN0 = deg2rad(RAAN0_deg);
aop0_deg = 50;
aop0 = deg2rad(aop0_deg);
ta0_deg = 0;
ta0 = deg2rad(ta0_deg);
tend = 30;

L = 15;%km
m1 = 250;%kg
m2 = 15;%kg
mt = 20;%kg
It = m2*(L*1000)^2 + mt*(1/3)*(L*1000)^2;
Ia = 0;
current_type = 2; % 0: constant, 2: controlled
current_val = 2;
tol = 1e-3;
dt = 100;

mag_model = 0; % [0: non-tilted; 1: igrf]

%% Inputs
menu1 = input("Run with defaults? [Y]/N", 's');
if isempty(menu1)
    menu1 = 'Y';
end
if menu1 == 'N'
    a0 = input("Semi-Major Axis [km] (default: "+num2str(a0)+"): ");
    if isempty(a0)
        a0 = 6878;
    end
    e0 = input("Eccentricity (default: "+num2str(e0)+"): ");
    if isempty(e0)
        e0 = .02;
    end
    i0_deg = input("Inclination [degree] (default: "+num2str(i0_deg)+"): ");
    if isempty(i0_deg)
        i0_deg = 30;
    end
    i0 = deg2rad(i0_deg);
    RAAN0_deg = input("Right Ascension of Ascending Node [degree] (default: "+num2str(RAAN0_deg)+"): ");
    if isempty(RAAN0_deg)
        RAAN0_deg = 30;
    end
    RAAN0 = deg2rad(RAAN0_deg);
    aop0_deg = input("Argument of Perigee [degree] (default: "+num2str(aop0_deg)+"): ");
    if isempty(aop0_deg)
        aop0_deg = 30;
    end
    aop0 = deg2rad(aop0_deg);
    ta0_deg = input("True Anomaly [degree] (default: "+num2str(ta0_deg)+"): ");
    if isempty(ta0_deg)
        ta0_deg = 30;
    end
    ta0 = deg2rad(ta0_deg);
    tend = input("Time [Hours] for simulation (default: "+num2str(tend)+"): ");
    if isempty(tend)
        tend = 30;
    end
    L = input("Tether length [km] (default "+num2str(L)+"): ");
    if isempty(L)
        L = 20;
    end
    m1 = input("Mass of main satellite [kg] (default: "+num2str(m1)+"): ");
    if isempty(m1)
        m1 = 250;
    end
    m2 = input("Mass of subsatellite [kg] (default: "+num2str(m2)+"): ");
    if isempty(m2)
        m2 = 150;
    end
    mt = input("Mass of tether [kg] (default: "+num2str(mt)+"): ");
    if isempty(mt)
        mt = 20;
    end
    It = input("Moment of inertia rotating the tether (default: "+num2str(It)+"): ");
    if isempty(It)
        It = m2*(L*1000)^2 + mt*(1/3)*(L*1000)^2;
    end
    Ia = input("Moment of inertia about tether axis (default: "+num2str(Ia)+"): ");
    if isempty(Ia)
        Ia = 0;
    end
    current_type = input("Type of Current control [0: constant; 1: OML; 2: controlled] (default: "+num2str(current_type)+"); ");
    if isempty(current_type)
        current_type = 0;
    end
    current_val = input("Value of Current [Value for constant; proportion for OML; max for control law] (default: "+num2str(current_val)+"); ");
    if isempty(current_val)
        current_val = 1;
    end
    mag_model = input("Magnetic model [0: non-tilted; 1: igrf] (default: "+num2str(mag_model)+"); ");
    if isempty(mag_model)
        mag_model = 0;
    end
    tol = input("Tolerance (dr/rp) [1e-2 is suggested by Vallado] (default: "+num2str(mag_model)+"); ");
    if isempty(mag_model)
        current_val = 1e-2;
    end
    dt = input("Time step [time is seconds] (default: "+num2str(mag_model)+"); ");
    if isempty(mag_model)
        current_val = 10;
    end
end
[ r0 , v0 ] = coes2state( [sqrt(396800*a0*(1-e0^2)), i0, e0, RAAN0, aop0, ta0], 396800 );
tspan = [0, tend*60*60];
sc_state0 = [r0; v0];
tether_state0 = [0.1; 0.1; 0; 0; 0; 0];
tether_param = [L; m1; m2; mt; It; It; Ia; current_type; current_val];

mu = 398600;
[ t , states] = encke_tether( tspan , sc_state0, tether_state0, ...
    tether_param, mag_model, mu , dt, tol );

%% Plot the Results
t = t/(60*60);
rvec = states(1:3,:);
vvec = states(4:6,:);
phi = rad2deg(states(7,:));
theta = rad2deg(states(8,:));
dphi = rad2deg(states(10,:));
dtheta = rad2deg(states(11,:));

lim_e = m2*L*1000*(1-cos(deg2rad(35)));
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

    I(ii) = current_val;
    pe_p(ii) = m2*L*1000*(1-cos(states(7,ii)));
    ke_p(ii) = (1/2)*It*(states(10,ii)*4)^2;
    if states(7,ii) < 0
        pe_p(ii) = -pe_p(ii);
    end
    if states(10,ii) < 0
        ke_p(ii) = -ke_p(ii);
    end
    e_p(ii) = pe_p(ii) + ke_p(ii);
    pe_t(ii) = m2*L*1000*(1-cos(states(8,ii)));
    ke_t(ii) = (1/2)*It*(states(11,ii)*4)^2;
    if states(8,ii) < 0
        pe_t(ii) = -pe_t(ii);
    end
    if states(11,ii) < 0
        ke_t(ii) = -ke_t(ii);
    end
    e_t(ii) = pe_t(ii) + ke_t(ii);
    if abs(e_p(ii)) > lim_e
        I(ii) = 0;
    elseif abs(e_t(ii)) > lim_e
        I(ii) = 0;
    end
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
subplot(2,1,1)
plot(t, theta)
xlabel('time [hours]')
ylabel('In Plane Libration [degree]')

subplot(2,1,2)
plot(t, phi)
xlabel('time [hours]')
ylabel('Out of Plane Libration [degree]')

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

figure
title("Current if Controlled")
subplot(3,1,1)
plot(t,I)
xlabel('time [hours]')
ylabel('Current [Amps]')

subplot(3,1,2)
plot(t,pe_p,t,ke_p,t,e_p,[t(1),tend(1)],[lim_e,lim_e],"r",[t(1),tend(1)],[-lim_e,-lim_e],"r")
xlabel('time [hours]')
ylabel('roll energy [J/m]')
legend("potential energy", "kinetic energy", "total energy","limit")

subplot(3,1,3)
plot(t,pe_t,t,ke_t,t,e_t,[t(1),tend(1)],[lim_e,lim_e],"r",[t(1),tend(1)],[-lim_e,-lim_e],"r")
xlabel('time [hours]')
ylabel('pitch energy [J/m]')
legend("potential energy", "kinetic energy", "total energy","limit")



