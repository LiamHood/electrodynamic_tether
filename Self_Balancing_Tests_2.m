clear; close all; clc;
% https://arc-aiaa-org.proxy.wichita.edu/doi/pdf/10.2514/6.2004-5309
% Self Balancing Tether
% m2 is upper mass
% m1 is lower mass
% mt is tether mass
% (m1, m2, mt) => (m, phi, LAMBDAt)
%% Define constants
m1 = 10; % lower mass
m2 = 10; % upper mass
mt = 10; % tether mass
L = 5000; % tether length
[m, phi, LAMBDAt, h_G, Is] = params_2_sbet_params(m1, m2, mt, L);
mue = 398600; % mu of earth
mue_si = mue*1e9; % mu of earth in meters
mum = 8e15; % magnetic constant of earths mag field

me = 9.1093837e-31; % mass of electron
mi = 2.6561e-26; % mass of atomic oxygen
e = 1.60217663e-19; % electron charge
mu = sqrt(me/mi); % Ratio (me/mi)^1/2
sigma = 3.5e7; % conductivity for aluminum

%% Define orbit
a = 7378; % [km], semi-major axis
e = .001; % eccentricity
i = deg2rad(25); % [degree], inclination
RAAN = deg2rad(0); % [degree], right ascension of ascending node
aop = deg2rad(0); % [degree], argument of perigee
ta = deg2rad(0); % [degree], true anomaly
coes_state0 = [a; e; i; RAAN; aop; ta];
[r, v] = coes2state([sqrt(mue*a*(1-e^2)),i,e,RAAN,aop,ta],mue);

%% Define attitude

% orbital frame is
% x is radial, y is momentum, z is backwards along track
% orbital frame definition in moving frame
ihato_m = r/norm(r);
jhato_m = cross(r,v)/norm(cross(r,v));
khato_m = cross(ihato_m,jhato_m);
orbital2moving = [ihato_m, jhato_m, khato_m];
moving2orbital = orbital2moving';

in_plane = deg2rad(10); % theta
out_plane = deg2rad(5); % phi
uhato = [cos(out_plane)*cos(in_plane); -sin(out_plane); cos(out_plane)*sin(in_plane)];
% [MG, ME] = sbet_torque(mue, a, i, ta, Is, orbital2moving*uhato, h_G, L);
tol = 1e-10;
tspan = [0, .2];
tether_state0 = [in_plane; out_plane; 0; 0];
opts = odeset('RelTol', tol, 'AbsTol', tol ) ;
[ t , states ] = ode45(@sbet_libration, tspan , [coes_state0; tether_state0], opts, mue) ;

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
xlabel("\phi, Out of Plane ")


function dstate = sbet_libration(t, state, mu)
    m1 = 10; % lower mass
    m2 = 10; % upper mass
    mt = 10; % tether mass
    L = 5000; % tether length
    [m, phi, LAMBDAt, h_G, Is] = params_2_sbet_params(m1, m2, mt, L);
    mue = 398600; % mu of earth
    mue_si = mue*1e9; % mu of earth in meters
    mum = 8e15; % magnetic constant of earths mag field
    
    % Tether Parameters
    m1 = 10; % lower mass
    m2 = 100; % upper mass
    mt = 10; % tether mass
    L = 5000; % tether length
    [m, phi, LAMBDAt, h_G, Is] = params_2_sbet_params(m1, m2, mt, L);
    tether_radius = .0005;
    p = tether_radius*2*pi; %perimeter in meters
    At = tether_radius^2*pi;
    
    % Constants
    me = 9.1093837e-31; % mass of electron
    mi = 2.6561e-26; % mass of atomic oxygen
    e = 1.60217663e-19; % electron charge
    mu = sqrt(me/mi); % Ratio (me/mi)^1/2
    sigma = 3.5e7; % conductivity for aluminum
    
    % Ionosphere/orbit parameters
    gamma = .15e-3; % secondary emission yield
    Em = .165; % Induced electric field; dot(cross(v, B), u)
    ninf = 2.0208e+11; % ionospheric plasma density 
    
    parameters = [p, At, me, Em, e, sigma, ninf, mu, gamma];
    
    % SOLVE BVP
    I_B_0 = 1;
    h0_0 = L/10;
    s0 = [0;1;I_B_0;h0_0];
    opts = optimoptions( 'fsolve' , 'Display' , 'off'  , 'FunctionTolerance' , 1e-8 ,...
        'StepTolerance' , 1e-8 , 'OptimalityTolerance' , 1e-8 , ...
        'Algorithm' , 'levenberg-marquardt') ; 
    [ s0s , F ] = fsolve( @(s)TetherSolveFun(s,parameters,L) , s0 , opts) ;
    h_span = [0, L];
    parameters = [p, At, me, Em, e, sigma, ninf, mu, gamma, s0s(4)];
    opts = odeset( 'RelTol' , 1e-8 , 'AbsTol' , 1e-8 ) ;
    [h, s] = ode45( @current_distribution_diffe, h_span, real(s0s(1:2)), opts, parameters) ;
    J1 = 0;
    for ii = 2:length(h)
        J1 = J1 + ((h_G-h(ii))*(s(ii,1))+(h_G-h(ii-1))*(s(ii-1,1)))/(h(ii)-h(ii-1));
    end
    epsilon = (J1/Is)*mum/mue_si;

    a = state(1);
    e = state(2);
    i = state(3);
    RAAN = state(4);
    aop = state(5);
    ta = state(6);
    theta = state(7); % in plane
    phi = state(8); % out of plane
    dtheta = state(9);
    dphi = state(10);

     % Calculate supplementary orbital elements
    p = a*(1-e^2);
    h = sqrt(mu*p);
    r = p/(1+e*cos(ta));
    n = sqrt(mu/a^3);
    u = aop + ta;

    f_roll = 0;
    f_pitch = 0;
    f_yaw = 0;
    da = (2/(n*sqrt(1-e^2)))*(e*sin(ta)*f_yaw + (p/r)*f_roll);
    de = (sqrt(1-e^2)/(n*a))*(sin(ta)*f_yaw + ...
        (cos(ta)+(e+cos(ta))/(1+e*cos(ta)))*f_roll);
    di = (r*cos(u)/(n*a^2*sqrt(1-e^2)))*f_pitch;
    dRAAN = (r*sin(u)/(n*a^2*sqrt(1-e^2)*sin(i)))*f_pitch;
    daop = (sqrt(1-e^2)/(n*a*e))*(-cos(ta)*f_yaw + ...
        sin(ta)*(1+r/p)*f_roll) - (r*cot(i)*sin(u)/h)*f_yaw;
    dta = h/r^2 + (1/(e*h))*(p*cos(ta)*f_yaw - (p+r)*sin(ta)*f_roll);
    ddtheta = 2*(1+dtheta)*dphi*tan(phi) - 1.5*sin(2*theta)...
        -epsilon*(sin(i)*tan(phi)*(2*sin(ta)*cos(theta) - cos(ta)*sin(theta)) + cos(i));
    ddphi = -sin(phi)*cos(phi)*((1+dtheta)^2 + 3*cos(theta)^2)...
        +epsilon*sin(i)*(2*sin(ta)*sin(theta)+cos(ta)*cos(theta));
    
    dstate = [da; de; di; dRAAN; daop; dta; dtheta; dphi; ddtheta; ddphi];
end