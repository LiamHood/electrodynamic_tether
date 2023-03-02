clear; close all; clc;

%% Defaults
mu = 398600; % [km^3/s^2]

% orbit
a0_def = 6878; % [km], semi-major axis
e0_def = .02; % eccentricity
i0_deg_def = 30; % [degree], inclination
RAAN0_deg_def = 30; % [degree], right ascension of ascending node
aop0_deg_def = 50; % [degree], argument of perigee
ta0_deg_def = 0; % [degree], true anomaly

% tether
L_def = 20; % [km], Tether length
m1_def = 750; % [kg], main satellite mass
m2_def = 250; % [kg], secondary satellite mass
mt_def = 20; % [kg], tether mass
Ia_def = 0; % [kg*m^2], inertia about local vertical axis
current_type_def = 1; % 0: constant, 2: controlled by energy limit
current_val_def = 5; % [A] current value, max if 


% simulation
tend_def = 30; % [hour], time of simulation
mag_model_def = 0; % 0: non-tilted; 1: igrf

%% Orbital inputs
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

%% Tether Inputs
menu2 = input("Adjust tether system? Y/[N]", 's');
fprintf("\tTether length [km]:"+num2str(L_def)+"\n")
fprintf("\tMass of main satellite [kg]: "+num2str(m1_def)+"\n");
fprintf("\tMass of subsatellite [kg]: "+num2str(m2_def)+"\n");
fprintf("\tMass of tether [kg]: "+num2str(mt_def)+"\n");
L_m = L_def; %tether length in meters
cm = (m1_def*0 + mt_def*L_m/2 + m2_def*L_m)/(m1_def+mt_def+m2_def); % find center of mass
I_m1 = m1_def*cm^2;
I_m2 = m2_def*(L_m-cm)^2;
I_mt = (1/12)*mt_def*L_m^2 + mt_def*(L_m/2 - cm)^2;
It_def = I_m1 + I_mt + I_m2; % kg*km^2
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
    It = input("Moment of inertia rotating the tether (Only input if different that barbell): ");
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
    L_m = L_def; %tether length in meters
    cm = (m1*0 + mt*L_m/2 + m2*L_m)/(m1+mt+m2); % find center of mass
    I_m1 = m1*cm^2;
    I_m2 = m2*(L_m-cm);
    I_mt = (1/12)*mt*L_m^2 + mt*(L_m/2 - cm);
    It = I_m1 + I_mt + I_m2;
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

%% Propagation Inputs
propagation_choice = input("Propagation method (V: variation of parameters, E: Encke's Method) [v]/E:","s");

if propagation_choice == "E" || propagation_choice == "e"
    % Encke's Method
    tol_def = 1e-4; % tolerance for propagation
    dt_def = 100; % [s], time step 

    % Inputs
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
    
    % Format inputs as needed for the simulation
    [ r0 , v0 ] = coes2state( [sqrt(mu*a0*(1-e0^2)), i0, e0, RAAN0, aop0, ta0], mu );
    tspan = [0, tend*60*60];
    sc_state0 = [r0; v0];
    tether_state0 = [0.0; 0.0; 0; 0; 0; 0];
    tether_param = [L; m1; m2; mt; It; It; Ia; current_type; current_val];
    
    % Run simulation using Encke's Method
    [ t , states] = encke_tether( tspan , sc_state0, tether_state0, ...
        tether_param, mag_model, mu , dt, tol );

    % Transform raw results for human readability
    t = t/(60*60);
    rvec = states(1:3,:);
    vvec = states(4:6,:);
    roll = rad2deg(states(7,:));
    pitch = rad2deg(states(8,:));
    droll = rad2deg(states(10,:));
    dpitch = rad2deg(states(11,:));
    
    % Determine result of energy controlled current and find COEs
    limit_libration = 35;
    I = current_val;
    m = m1+m2+mt;
    lim_e = m*L*(1-cos(deg2rad(limit_libration)));
    for ii = 1:length(t)
        COES = state2coes_display([rvec(:,ii), vvec(:,ii)], mu );
        a(ii) = COES(9);
        e(ii) = COES(3);
        i(ii) = COES(2);
        RAAN(ii) = COES(4);
        aop(ii) = COES(5);
        ta(ii) = COES(6);
        r(ii) = norm(rvec(:,ii));
        v(ii) = norm(vvec(:,ii));
    
        I(ii) = current_val;
        I(ii) = current_val;
        pe_roll(ii) = m*L*(1-cos(roll(ii)));
        ke_roll(ii) = (1/2)*m*(L*droll(ii))^2;
        if roll(ii) < 0
            pe_roll(ii) = -pe_roll(ii);
        end
        if droll(ii) < 0
            ke_roll(ii) = -ke_roll(ii);
        end
        e_roll(ii) = pe_roll(ii) + ke_roll(ii);
        pe_pitch(ii) = m*L*(1-cos(pitch(ii)));
        ke_pitch(ii) = (1/2)*m*(L*dpitch(ii))^2;
%             ke_t = (1/2)*Iy*(dtheta)^2;
        if pitch(ii) < 0
            pe_pitch(ii) = -pe_pitch(ii);
        end
        if dpitch(ii) < 0
            ke_pitch(ii) = -ke_pitch(ii);
        end
        e_pitch(ii) = pe_pitch(ii) + ke_pitch(ii);
        if abs(e_roll(ii)) > lim_e
            I(ii) = 0;
        elseif abs(e_pitch(ii)) > lim_e
            I(ii) = 0;
        end
    end
else
    % Variation of Parameters
    tol_def = 1e-8;

    % Inputs
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

    % Format inputs as needed for the simulation
    tspan = [0, tend*60*60];
    sc_state0 = [a0; e0; i0; RAAN0; aop0; ta0];
    tether_state0 = [0.1; 0.1; 0; 0; 0; 0];
    tether_param = [L; m1; m2; mt; It; It; Ia; current_type; current_val; .0006];
    
    % Run simulation using Variation of Parameters
    [ t , states] = BasicTether( tspan , sc_state0, tether_state0, ...
        tether_param, mu , tol );
    
    % Transform raw results for human readability
    t = t/(60*60);
    a = states(:,1);
    e = states(:,2);
    i = rad2deg(states(:,3));
    RAAN = rad2deg(states(:,4));
    aop = rad2deg(states(:,5));
    ta = rad2deg(states(:,6));
    
    roll = rad2deg(states(:,7));
    pitch = rad2deg(states(:,8));
    yaw = rad2deg(states(:,9));
    droll = rad2deg(states(:,10));
    dpitch = rad2deg(states(:,11));
    
    % Determine result of energy controlled current
    limit_libration = 35;
    I = current_val;
    m = m1+m2+mt;
    lim_e = L*1000*(1-cos(deg2rad(limit_libration)));
%     lim_e = limit_libration;
    for ii = 1:length(t)
        [ rvec(:,ii) , vvec(:,ii) ] = coes2state( [sqrt(mu*a(ii)*(1-e(ii)^2)), ...
            states(ii,3), e(ii), states(ii,4), states(ii,5), states(ii,6)] , mu );
    
        r(ii) = norm(rvec(:,ii));
        v(ii) = norm(vvec(:,ii));
    
        [I_c, db] = controlled_current(limit_libration, current_val, states(ii,7), states(ii,8), states(ii,10), states(ii,11), m, L*1000);
        I(ii) = I_c;
        pe_roll(ii) = db(2);
        ke_roll(ii) = db(3);
        e_roll(ii) = db(4);
        pe_pitch(ii) = db(5);
        ke_pitch(ii) = db(6);
        e_pitch(ii) = db(7);

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
plot(t, pitch, t, dpitch.^2*6e4, t, abs(pitch)+abs(dpitch.^2)*6e4)
xlabel('time [hours]')
ylabel('In Plane Libration [degree]')
legend('displacement', 'rate', 'energy')

subplot(2,1,2)
plot(t, roll, t, droll.^2*6e4, t, abs(roll)+abs(droll.^2)*6e4)
xlabel('time [hours]')
ylabel('Out of Plane Libration [degree]')

figure
subplot(2,1,1)
plot(t, pitch)
xlabel('time [hours]')
ylabel('In Plane Libration [degree]')
legend('displacement', 'rate', 'energy')

subplot(2,1,2)
plot(t, roll)
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
plot(t,pe_roll,t,ke_roll,t,e_roll,[t(1),tend(1)],[lim_e,lim_e],"r",[t(1),tend(1)],[-lim_e,-lim_e],"r")
xlabel('time [hours]')
ylabel('roll energy [J/m]')
legend("potential energy", "kinetic energy", "total energy","limit")

subplot(3,1,3)
plot(t,pe_pitch,t,ke_pitch,t,e_pitch,[t(1),tend(1)],[lim_e,lim_e],"r",[t(1),tend(1)],[-lim_e,-lim_e],"r")
xlabel('time [hours]')
ylabel('pitch energy [J/m]')
legend("potential energy", "kinetic energy", "total energy","limit")



