clear; close all; clc;
% https://arc-aiaa-org.proxy.wichita.edu/doi/pdf/10.2514/6.2004-5309
% Self Balancing Tether
% m2 is upper mass
% m1 is lower mass
% mt is tether mass
% (m1, m2, mt) => (m, phi, LAMBDAt)
%% Define tether
tether.m1 = 321.6; % lower mass [kg]
tether.m2 = 563.5; % upper mass [kg]
tether.L = 20000; % tether length [km]
tether.density = 2700; % [kg/m^3]
tether.thickness = .18e-3; % [m]
tether.perimeter = 24e-3; % [m]
tether.At = (tether.perimeter-tether.thickness*2)/2*tether.thickness; % [m^2]
tether.mt = tether.density*tether.L*tether.At; % tether mass [kg]
[tether.m, tether.phi, tether.LAMBDAt, tether.h_G, tether.Is] = params_2_sbet_params(tether.m1, tether.m2, tether.mt, tether.L); % [kg, radians, non-dimensional, m, kg*m^2]
tether.sigma = 3.5e7; % conductivity for aluminum [1/(Ohm*m)]

%% define constants
constants.mue = 398600; % mu of earth [km^3/s^2]
constants.mue_si = constants.mue*1e9; % mu of earth in meters [m^3/s^2]
A_m2_mum = 8.22e22; % magnetic constant of earths mag field [A*m^2] https://www.britannica.com/science/geomagnetic-field/Dipolar-field
constants.mum = A_m2_mum/795774.71545947673925; % magnetic constant of earths mag field [T*M^3] https://maurermagnetic.com/en/demagnetizing/technology/convert-magnetic-units/
constants.me = 9.1093837e-31; % mass of electron [kg]
constants.mi = 2.6561e-26; % mass of atomic oxygen [kg]
constants.e = 1.60217663e-19; % electron charge [C]
constants.mu = sqrt(constants.me/constants.mi); % Ratio 
constants.gamma = .15e-3; % secondary emission yield [1/V]

%% Define orbit
a = 7378; % [km], semi-major axis
ecc = .001; % eccentricity
i = deg2rad(25); % [radian], inclination
RAAN = deg2rad(0); % [radian], right ascension of ascending node
aop = deg2rad(0); % [radian], argument of perigee
ta = deg2rad(0); % [radian], true anomaly
coes_state0 = [a; ecc; i; RAAN; aop; ta];
[r, v] = coes2state([sqrt(constants.mue*a*(1-ecc^2)),i,ecc,RAAN,aop,ta],constants.mue);

%% Define attitude
% C_state2orbital = state_2_SBETorbital(r, v);

in_plane = deg2rad(0); % theta [radian]
out_plane = deg2rad(0); % phi [radian]



tol = 1e-6;
tspan = [0, 10*3600];
tether_state0 = [in_plane; out_plane; 0; 0];
opts = odeset('RelTol', tol, 'AbsTol', tol ) ;
e_density = load("nna_llat2densityL5N8.mat", "net");
bvp_ann = load("nna_L1N128_1e2_1e2.mat");
[ t , states ] = ode45(@sbet_libration_ann, tspan , [coes_state0; tether_state0], opts, constants, tether, e_density.net, bvp_ann.net) ;


% Find EM
for ii = 1:length(t)
    state = states(ii,:);
    theta = states(ii,7);
    phi = states(ii,8);

    a = state(1);
    ecc = state(2);
    i = state(3);
    RAAN = state(4);
    aop = state(5);
    ta = state(6);
    theta = state(7); % in plane
    phi = state(8); % out of plane
    dtheta = state(9);
    dphi = state(10);

     % Calculate supplementary orbital elements
    p = a*(1-ecc^2);
    h = sqrt(constants.mue*p);
    rmag = p/(1+ecc*cos(ta));
    n = sqrt(constants.mue/a^3);
    u = aop + ta;
    [ r , v ] = coes2state([h, i, ecc, RAAN, aop, ta], constants.mue);

    [Bx, By, Bz] = MagField_OrbitalFrame(state, constants.mum); % a, e, i, RAAN, omega, theta
    B = [Bx; By; Bz];
    tether_rotation = tether_attitude(theta, phi);
    uhat = tether_rotation*[1;0;0];
    C_state2orbital = state_2_SBETorbital(r, v);
    v_orbit = C_state2orbital*v*1e3;
    Em(ii) = dot(cross(v_orbit, B), uhat);
end

figure
hold on
plot(t, rad2deg(states(:, 7)))
plot(t, rad2deg(states(:, 8)))
hold off
legend("\theta, In Plane", "\phi, Out of Plane")
xlabel("Time [s]")
ylabel("Libration ")

figure
plot(rad2deg(states(:, 7)), rad2deg(states(:, 8)))
xlabel("\theta, In Plane ")
ylabel("\phi, Out of Plane ")

figure
plot(t, states(:,1))
xlabel("Time [s]")
ylabel("Semi-major Axis [km]")

figure
plot(t, Em)
xlabel("Time [s]")
ylabel("Induced Electric Field [V/m]")

function dstate = sbet_libration(t, state, constants, tether, e_density_ann)
    t/3600
    jdate = 2458849.5 + t/(24*3600);

    a = state(1);
    ecc = state(2);
    i = state(3);
    RAAN = state(4);
    aop = state(5);
    ta = state(6);
    theta = state(7); % in plane
    phi = state(8); % out of plane
    dtheta = state(9);
    dphi = state(10);

     % Calculate supplementary orbital elements
    p = a*(1-ecc^2);
    h = sqrt(constants.mue*p);
    rmag = p/(1+ecc*cos(ta));
    n = sqrt(constants.mue/a^3);
    u = aop + ta;
    [ r , v ] = coes2state([h, i, ecc, RAAN, aop, ta], constants.mue);

    parameters.ninf = find_electron_density(jdate, r, e_density_ann);

    [Bx, By, Bz] = MagField_OrbitalFrame(state, constants.mum); % a, e, i, RAAN, omega, theta
    B = [Bx; By; Bz];
    tether_rotation = tether_attitude(theta, phi);
    uhat = tether_rotation*[1;0;0];
    C_state2orbital = state_2_SBETorbital(r, v);
    v_orbit = C_state2orbital*v*1e3;
    parameters.Em = dot(cross(v_orbit, B), uhat);
    
    % SOLVE BVP
    h0_0 = tether.L/10;
    s0 = [1;h0_0];
    opts = optimoptions( 'fsolve' , 'Display' , 'none'  , 'FunctionTolerance' , 1e-8 ,...
        'StepTolerance' , 1e-8 , 'OptimalityTolerance' , 1e-8 , ...
        'Algorithm' , 'levenberg-marquardt') ; 
    [ s0s , F ] = fsolve( @(s)TetherSolveFun2(s, tether, constants, parameters) , s0 , opts) ;
    PHIe_0 = s0s(1); 
    h_span = [0, tether.L];
    h_0 = s0s(2);
    opts = odeset( 'RelTol' , 1e-10 , 'AbsTol' , 1e-10 ) ;
    [h, s] = ode45( @current_distribution_diffe2, h_span, [0, PHIe_0], opts, h_0, tether, constants, parameters) ;
    J1 = 0;
    Iint = 0;
    fe_right = 0;
    for ii = 2:length(h)
        J1 = J1 + ((tether.h_G-h(ii))*(s(ii,1))+(tether.h_G-h(ii-1))*(s(ii-1,1)))/2*(h(ii)-h(ii-1));
        Iint = Iint + (s(ii,1)+s(ii-1,1))/2*(h(ii)-h(ii-1));
        fe_right = fe_right + (h(ii)*s(ii,1) + h(ii-1)*s(ii-1,1))/2 * (h(ii)-h(ii-1));
    end
    epsilon = (J1/tether.Is)*constants.mum/constants.mue_si;
    Im = Iint/tether.L;
    fe = fe_right/(tether.L^2*Im);

    % calculate force
    f_orbital = Im*tether.L*cross(uhat, [Bx, By, Bz]);
    a_orbital = f_orbital/tether.m*1e-3'; % force to acceleraction in m to km
    a_r = a_orbital(1);
    a_h = a_orbital(2);
    a_theta = a_orbital(3);

    hmag = sqrt(constants.mue*p);
    da = (2/(n*sqrt(1-ecc^2)))*(ecc*sin(ta)*a_r + (p/rmag)*a_theta);
    de = (sqrt(1-ecc^2)/(n*a))*(sin(ta)*a_r + ...
        (cos(ta)+(ecc+cos(ta))/(1+ecc*cos(ta)))*a_theta);
    di = (rmag*cos(u)/(n*a^2*sqrt(1-ecc^2)))*a_h;
    dRAAN = (rmag*sin(u)/(n*a^2*sqrt(1-ecc^2)*sin(i)))*a_h;
    daop = (sqrt(1-ecc^2)/(n*a*ecc))*(-cos(ta)*a_r + ...
        sin(ta)*(1+rmag/p)*a_theta) - (rmag*cot(i)*sin(u)/hmag)*a_r;
    dta = hmag/rmag^2 + (1/(ecc*hmag))*(p*cos(ta)*a_r - (p+rmag)*sin(ta)*a_theta);
    
    Tq = edt_torque_orbital(theta, phi, tether.L, fe, Im, Bx, By, Bz);
    Tgg = gravity_grad_torque_orbital(rmag*1e3, theta, phi, tether.Is, constants.mue_si);
    T = Tq+Tgg;
    % x is towards zenith; 0
    % y is -h, rotation is in-plane or theta; Is
    % z is along track, rotation is out-of-plane or phi; Is
    ddtheta = T(2)/tether.Is;
    ddphi = T(3)/tether.Is;

    dstate = [da; de; di; dRAAN; daop; dta; dtheta; dphi; ddtheta; ddphi];
%     disp(dstate')
end

function dstate = sbet_libration_ann(t, state, constants, tether, e_density_ann, bvp_ann)
    t/3600
    jdate = 2458849.5 + t/(24*3600);

    a = state(1);
    ecc = state(2);
    i = state(3);
    RAAN = state(4);
    aop = state(5);
    ta = state(6);
    theta = state(7); % in plane
    phi = state(8); % out of plane
    dtheta = state(9);
    dphi = state(10);

     % Calculate supplementary orbital elements
    p = a*(1-ecc^2);
    h = sqrt(constants.mue*p);
    rmag = p/(1+ecc*cos(ta));
    n = sqrt(constants.mue/a^3);
    u = aop + ta;
    [ r , v ] = coes2state([h, i, ecc, RAAN, aop, ta], constants.mue);

    parameters.ninf = find_electron_density(jdate, r, e_density_ann);

    [Bx, By, Bz] = MagField_OrbitalFrame(state, constants.mum); % a, e, i, RAAN, omega, theta
    B = [Bx; By; Bz];
    tether_rotation = tether_attitude(theta, phi);
    uhat = tether_rotation*[1;0;0];
    C_state2orbital = state_2_SBETorbital(r, v);
    v_orbit = C_state2orbital*v*1e3;
    parameters.Em = dot(cross(v_orbit, B), uhat);
    
    % SOLVE BVP
%     h0_0 = tether.L/10;
%     s0 = [1;h0_0];
%     opts = optimoptions( 'fsolve' , 'Display' , 'none'  , 'FunctionTolerance' , 1e-8 ,...
%         'StepTolerance' , 1e-8 , 'OptimalityTolerance' , 1e-8 , ...
%         'Algorithm' , 'levenberg-marquardt') ; 
%     [ s0s , F ] = fsolve( @(s)TetherSolveFun2(s, tether, constants, parameters) , s0 , opts) ;
%     PHIe_0 = s0s(1); 
%     h_span = [0, tether.L];
%     h_0 = s0s(2);

    % Use ANN to approximate BVP
    ann_input = [1; 1; 1; parameters.Em/4; parameters.ninf/2e12];
    ann_out = bvp_ann(ann_input);
    PHIe_0 = ann_out(1)*1e5;
    h_0 = ann_out(2)*1e3;
    h_span = [0, tether.L];
    opts = odeset( 'RelTol' , 1e-10 , 'AbsTol' , 1e-10 ) ;
    [h, s] = ode45( @current_distribution_diffe2, h_span, [0, PHIe_0], opts, h_0, tether, constants, parameters) ;
    J1 = 0;
    Iint = 0;
    fe_right = 0;
    for ii = 2:length(h)
        J1 = J1 + ((tether.h_G-h(ii))*(s(ii,1))+(tether.h_G-h(ii-1))*(s(ii-1,1)))/2*(h(ii)-h(ii-1));
        Iint = Iint + (s(ii,1)+s(ii-1,1))/2*(h(ii)-h(ii-1));
        fe_right = fe_right + (h(ii)*s(ii,1) + h(ii-1)*s(ii-1,1))/2 * (h(ii)-h(ii-1));
    end
    epsilon = (J1/tether.Is)*constants.mum/constants.mue_si;
    Im = Iint/tether.L;
    fe = fe_right/(tether.L^2*Im);

    % calculate force
    f_orbital = Im*tether.L*cross(uhat, [Bx, By, Bz]);
    a_orbital = f_orbital/tether.m*1e-3'; % force to acceleraction in m to km
    a_r = a_orbital(1);
    a_h = a_orbital(2);
    a_theta = a_orbital(3);

    hmag = sqrt(constants.mue*p);
    da = (2/(n*sqrt(1-ecc^2)))*(ecc*sin(ta)*a_r + (p/rmag)*a_theta);
    de = (sqrt(1-ecc^2)/(n*a))*(sin(ta)*a_r + ...
        (cos(ta)+(ecc+cos(ta))/(1+ecc*cos(ta)))*a_theta);
    di = (rmag*cos(u)/(n*a^2*sqrt(1-ecc^2)))*a_h;
    dRAAN = (rmag*sin(u)/(n*a^2*sqrt(1-ecc^2)*sin(i)))*a_h;
    daop = (sqrt(1-ecc^2)/(n*a*ecc))*(-cos(ta)*a_r + ...
        sin(ta)*(1+rmag/p)*a_theta) - (rmag*cot(i)*sin(u)/hmag)*a_r;
    dta = hmag/rmag^2 + (1/(ecc*hmag))*(p*cos(ta)*a_r - (p+rmag)*sin(ta)*a_theta);
    
    Tq = edt_torque_orbital(theta, phi, tether.L, fe, Im, Bx, By, Bz);
    Tgg = gravity_grad_torque_orbital(rmag*1e3, theta, phi, tether.Is, constants.mue_si);
    T = Tq+Tgg;
    % x is towards zenith; 0
    % y is -h, rotation is in-plane or theta; Is
    % z is along track, rotation is out-of-plane or phi; Is
    ddtheta = T(2)/tether.Is;
    ddphi = T(3)/tether.Is;

    dstate = [da; de; di; dRAAN; daop; dta; dtheta; dphi; ddtheta; ddphi];
%     disp(dstate')
end
