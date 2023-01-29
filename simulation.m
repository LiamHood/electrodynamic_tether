clear; close all; clc;

%% Defaults
mu = 398600;

% orbit
a0_def = 6978; % [km], semi-major axis
e0_def = .02; % eccentricity
i0_deg_def = 30; % [degree], inclination
RAAN0_deg_def = 30; % [degree], right ascension of ascending node
aop0_deg_def = 50; % [degree], argument of perigee
ta0_deg_def = 0; % [degree], true anomaly

% tether
L_def = 15; % [km], Tether length
m1_def = 250; % [kg], main satellite mass
m2_def = 15; % [kg], secondary satellite mass
mt_def = 20; % [kg], tether mass
Ia_def = 0; % [kg*m^2], inertia about local vertical axis
current_type_def = 0; % 0: constant, 2: controlled by energy limit
current_val_def = 1; % [A] current value, max if 


% simulation
tend_def = 30;
tol_def = 1e-4;
dt_def = 100;
mag_model_def = 0; % [0: non-tilted; 1: igrf]

%% Inputs
menu1 = input("Adjust orbit initial conditions? Y/[N]", 's');
fprintf("\tSemi-Major Axis [km]: "+num2str(a0_def)+"\n")
fprintf("\tEccentricity: "+num2str(e0_def)+"\n");
fprintf("\tInclination [degree]: "+num2str(i0_deg_def)+"\n");
fprintf("\tRight Ascension of Ascending Node [degree]: "+num2str(RAAN0_deg_def)+"\n");
fprintf("\tArgument of Perigee [degree]: "+num2str(aop0_deg_def)+"\n");
fprintf("\tTrue Anomaly [degree]: "+num2str(ta0_deg_def)+"\n");
if isempty(menu1)
    menu1 = 'N';
end
if menu1 == 'Y'
    a0 = input("Semi-Major Axis [km]: ");
    e0 = input("Eccentricity: ");
    i0_deg = input("Inclination [degree]: ");
    RAAN0_deg = input("Right Ascension of Ascending Node [degree]: ");
    aop0_deg = input("Argument of Perigee [degree]: ");
    ta0_deg = input("True Anomaly [degree]: ");
else 
    a0 = [];
    e0 = [];
    i0_deg = [];
    RAAN0_deg = [];
    aop0_deg = [];
    ta0_deg = [];
end
if isempty(a0)
    a0 = a0_def;
end
if isempty(e0)
    e0 = e0_def;
end
if isempty(i0_deg)
    i0_deg = i0_deg_def;
end
i0 = deg2rad(i0_deg);
if isempty(RAAN0_deg)
    RAAN0_deg = RAAN0_deg_def;
end
RAAN0 = deg2rad(RAAN0_deg);
if isempty(aop0_deg)
    aop0_deg = aop0_deg_def;
end
aop0 = deg2rad(aop0_deg);
if isempty(ta0_deg)
    ta0_deg = ta0_deg_def;
end
ta0 = deg2rad(ta0_deg);


menu2 = input("Adjust tether system? Y/[N]", 's');
fprintf("\tTether length [km]:"+num2str(L_def)+"\n")
fprintf("\tMass of main satellite [kg]: "+num2str(m1_def)+"\n");
fprintf("\tMass of subsatellite [kg]: "+num2str(m2_def)+"\n");
fprintf("\tMass of tether [kg]: "+num2str(mt_def)+"\n");
It_def = m2_def*(L_def*1000)^2 + mt_def*(1/3)*(L_def*1000)^2;
fprintf("\tMoment of inertia rotating the tether: "+num2str(It_def)+"\n");
fprintf("\tMoment of inertia about tether axis: "+num2str(Ia_def)+"\n");
fprintf("\tType of Current control [0: constant; 1: OML limited; 2: controlled]: "+num2str(current_type_def)+"\n");
fprintf("\tValue of Current [Value for constant; proportion for OML; max for controlled]: "+num2str(current_val_def)+"\n");
if isempty(menu2)
    menu1 = 'N';
end
if menu2 == 'Y'
    L = input("Tether length [km]: ");
    m1 = input("Mass of main satellite [kg]: ");
    m2 = input("Mass of subsatellite [kg]: ");
    mt = input("Mass of tether [kg]: ");
    It = input("Moment of inertia rotating the tether: ");
    Ia = input("Moment of inertia about tether axis: ");
    current_type = input("Type of Current control [0: constant; 1: OML limited; 2: controlled]: ");
    current_val = input("Value of Current [Value for constant; proportion for OML; max for control law]: ");
else
    L = [];
    m1 = [];
    m2 = [];
    mt = [];
    It = [];
    Ia = [];
    current_type = [];
    current_val = [];
end
 if isempty(L)
    L = L_def;
 end
 if isempty(m1)
    m1 = m1_def;
end
if isempty(m2)
    m2 = m2_def;
end
if isempty(mt)
    mt = mt_def;
end
if isempty(It)
    It = m2*(L*1000)^2 + mt*(1/3)*(L*1000)^2;
end
if isempty(Ia)
    Ia = Ia_def;
end
if isempty(current_type)
    current_type = current_type_def;
end
if isempty(current_val)
    current_val = current_val_def;
end


propagation_choice = input("Propagation method (V: variation of parameters, E: Encke's Method) [v]/E:","s");

if propagation_choice == "E" || propagation_choice == "e"
    menu3 = input("Adjust simulation settings? Y/[N]", 's');
    fprintf("\tMagnetic model [0: non-tilted; 1: igrf]: "+num2str(mag_model_def)+"\n")
    fprintf("\tTolerance (dr/rp) [1e-2 is suggested by Vallado]: "+num2str(mag_model_def)+"\n");
    fprintf("\tTime step [time is seconds]: "+num2str(mag_model_def)+"\n");
    fprintf("\tTime [Hours] for simulation: "+num2str(tend_def)+"\n");
    if isempty(menu3)
        menu1 = 'N';
    end
    if menu3 == 'Y'
        mag_model = input("Magnetic model [0: non-tilted; 1: igrf]: ");
        tol = input("Tolerance (dr/rp): ");
        dt = input("Time step [time is seconds]: ");
        tend = input("Time [Hours] for simulation: ");
    else
        mag_model = [];
        tol = [];
        dt = [];
        tend = [];
    end
    if isempty(mag_model)
        mag_model = mag_model_def;
    end
    if isempty(tol)
        tol = tol_def;
    end
    if isempty(dt)
        dt = dt_def;
    end
    if isempty(tend)
        tend = tend_def;
    end
    
    [ r0 , v0 ] = coes2state( [sqrt(mu*a0*(1-e0^2)), i0, e0, RAAN0, aop0, ta0], mu );
    tspan = [0, tend*60*60];
    sc_state0 = [r0; v0];
    tether_state0 = [0.1; 0.1; 0; 0; 0; 0];
    tether_param = [L; m1; m2; mt; It; It; Ia; current_type; current_val];
    
    [ t , states] = encke_tether( tspan , sc_state0, tether_state0, ...
        tether_param, mag_model, mu , dt, tol );

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
else
    tol_def = 1e-8;
    menu3 = input("Adjust simulation settings? Y/[N]", 's');
    fprintf("\tMagnetic model [0: non-tilted; 1: igrf]: "+num2str(mag_model_def)+"\n")
    fprintf("\tTolerance [for ode45]: "+num2str(mag_model_def)+"\n");
    fprintf("\tTime [Hours] for simulation: "+num2str(tend_def)+"\n");
    if isempty(menu3)
        menu1 = 'N';
    end
    if menu3 == 'Y'
        mag_model = input("Magnetic model [0: non-tilted; 1: igrf]: ");
        tol = input("Tolerance [for ode45]: ");
        tend = input("Time [Hours] for simulation: ");
    else
        mag_model = [];
        tol = [];
        dt = [];
        tend = [];
    end
    if isempty(mag_model)
        mag_model = mag_model_def;
    end
    if isempty(tol)
        tol = tol_def;
    end
    if isempty(tend)
        tend = tend_def;
    end

    tspan = [0, tend*60*60];
    sc_state0 = [a0; e0; i0; RAAN0; aop0; ta0];
    tether_state0 = [0.1; 0.1; 0; 0; 0; 0];
    tether_param = [L; m1; m2; mt; It; It; Ia; current_type; current_val];
    
    [ t , states] = BasicTether( tspan , sc_state0, tether_state0, ...
        tether_param, mu , tol );
    
    t = t/(60*60);
    a = states(:,1);
    e = states(:,2);
    i = rad2deg(states(:,3));
    RAAN = rad2deg(states(:,4));
    aop = rad2deg(states(:,5));
    ta = rad2deg(states(:,6));
    
    phi = rad2deg(states(:,7));
    theta = rad2deg(states(:,8));
    psi = rad2deg(states(:,9));
    dphi = rad2deg(states(:,10));
    dtheta = rad2deg(states(:,11));
    
    lim_e = m2*L*1000*(1-cos(deg2rad(35)));
    for ii = 1:length(t)
        [ rvec(:,ii) , vvec(:,ii) ] = coes2state( [sqrt(mu*a(ii)*(1-e(ii)^2)), ...
            states(ii,3), e(ii), states(ii,4), states(ii,5), states(ii,6)] , mu );
    
        r(ii) = norm(rvec(:,ii));
        v(ii) = norm(vvec(:,ii));
    
        I(ii) = current_val;
        pe_p(ii) = m2*L*1000*(1-cos(states(ii,7)));
        ke_p(ii) = (1/2)*It*(states(ii,10)*4)^2;
        if states(ii,7) < 0
            pe_p(ii) = -pe_p(ii);
        end
        if states(ii,10) < 0
            ke_p(ii) = -ke_p(ii);
        end
        e_p(ii) = pe_p(ii) + ke_p(ii);
        pe_t(ii) = m2*L*1000*(1-cos(states(ii,8)));
        ke_t(ii) = (1/2)*It*(states(ii,11)*4)^2;
        if states(ii,8) < 0
            pe_t(ii) = -pe_t(ii);
        end
        if states(ii,11) < 0
            ke_t(ii) = -ke_t(ii);
        end
        e_t(ii) = pe_t(ii) + ke_t(ii);
        if abs(e_p(ii)) > lim_e
            I(ii) = 0;
        elseif abs(e_t(ii)) > lim_e
            I(ii) = 0;
        end
    
    end
end

%% Plot the Results

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



