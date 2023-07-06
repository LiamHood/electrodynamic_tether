clear; close all; clc;
% https://arc-aiaa-org.proxy.wichita.edu/doi/pdf/10.2514/6.2004-5309
% Self Balancing Tether
% m2 is upper mass
% m1 is lower mass
% mt is tether mass
% (m1, m2, mt) => (m, phi, LAMBDAt)
%% Define constants
m1 = 100; % lower mass
m2 = 200; % upper mass
mt = 50; % tether mass
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
jhato_m = -cross(r,v)/norm(cross(r,v));
khato_m = cross(ihato_m,jhato_m)/norm(cross(ihato_m,jhato_m));
orbital2moving = [ihato_m, jhato_m, khato_m];
moving2orbital = orbital2moving';

in_plane = deg2rad(0); % theta
out_plane = deg2rad(0); % phi
u1 = [cos(out_plane)*cos(in_plane); -sin(out_plane); cos(out_plane)*sin(in_plane)];
u2 = [-sin(in_plane); 0; cos(in_plane)];
u3 = [-sin(out_plane)*cos(in_plane); -cos(out_plane); -sin(out_plane)*sin(in_plane)];
tether_rotation = [u1, u2, u3];
uhato = tether_rotation*[1;0;0];
tol = 1e-6;
tspan = [0, 6];
tether_state0 = [in_plane; out_plane; 0; 0];
opts = odeset('RelTol', tol, 'AbsTol', tol ) ;
load("nna_llat2densityL5N6.mat")
[ t , states ] = ode45(@sbet_libration, tspan , [coes_state0; tether_state0], opts, mue, net) ;

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


function dstate = sbet_libration(t, state, mue, net)
    t
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
    h = sqrt(mue*p);
    r = p/(1+ecc*cos(ta));
    n = sqrt(mue/a^3);
    u = aop + ta;

    % Tether Parameters
    m1 = 100; % lower mass
    m2 = 200; % upper mass
    mt = 50; % tether mass
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
    mum = 8e15; % magnetic constant of earths mag field
    mue_si = mue*1e9; % mu of earth in meters

    % Ionosphere/orbit parameters
    gamma = .15e-3; % secondary emission yield
    
%     Em = .165; % Induced electric field; dot(cross(v, B), u)
    
    utc_time = datetime(jdate,'ConvertFrom','juliandate');
    [ r , v ] = coes2state([sqrt(mue*a*(1-e^2)), i, e, RAAN, aop, ta], mue );
    lla = eci2lla(r'*1e3,[year(utc_time), month(utc_time), day(utc_time),...
        hour(utc_time), minute(utc_time),second(utc_time)]);
    UT = hour(utc_time)+minute(utc_time)/60+second(utc_time)/3600;
    if lla(2)>0
        LT = UT + lla(2)*(12/180);
    else
        LT = UT + (lla(2)+360)*(12/180);
    end
    LT_norm = LT/24;
    lat_norm = lla(1)/90;
    long_norm = lla(2)/180;
    alt_norm = (lla(3)/1e3)/2000; % m to km, then normalize to neural net
    ninf = net([lat_norm; long_norm; alt_norm; LT_norm])*2e12; % ionospheric plasma density 
    [Bx, By, Bz] = MagField_OrbitalFrame(state); % a, e, i, RAAN, omega, theta
    B = [Bx; By; Bz];
    uhat = [cos(phi)*cos(theta), -sin(phi), cos(phi)*sin(theta)]; % paper config
    i_orbit = r./norm(r);
    j_orbit = -cross(r,v)./norm(cross(r,v));
    k_orbit = cross(i_orbit, j_orbit)./norm(cross(i_orbit, j_orbit));
    orbit2eci = [i_orbit, j_orbit, k_orbit];
    v_orbit = orbit2eci'*v*1e3;
    Em = dot(cross(v_orbit, B), uhat);
    parameters = [p, At, me, Em, e, sigma, ninf, mu, gamma];
    
    % SOLVE BVP
    I_B_0 = 1;
    h0_0 = L/10;
    s0 = [1;I_B_0;h0_0];
    opts = optimoptions( 'fsolve' , 'Display' , 'off'  , 'FunctionTolerance' , 1e-8 ,...
        'StepTolerance' , 1e-8 , 'OptimalityTolerance' , 1e-8 , ...
        'Algorithm' , 'levenberg-marquardt') ; 
    [ s0s , F ] = fsolve( @(s)TetherSolveFun(s,parameters,L) , s0 , opts) ;
    PHIe_0 = s0s(1); 
    I_B = s0s(2);
    h_span = [0, L];
    h_0 = s0s(3);
    parameters = [p, At, me, Em, e, sigma, ninf, mu, gamma, h_0];
    opts = odeset( 'RelTol' , 1e-8 , 'AbsTol' , 1e-8 ) ;
    [h, s] = ode45( @current_distribution_diffe, h_span, [0, PHIe_0], opts, parameters) ;
    J1 = 0;
    Iint = 0;
    fe_right = 0;
    for ii = 2:length(h)
        J1 = J1 + ((h_G-h(ii))*(s(ii,1))+(h_G-h(ii-1))*(s(ii-1,1)))/2*(h(ii)-h(ii-1));
        Iint = Iint + (s(ii,1)+s(ii-1,1))/2*(h(ii)-h(ii-1));
        fe_right = fe_right + (h(ii)*s(ii,1) + h(ii-1)*s(ii-1,1))/2 * (h(ii)-h(ii-1));
    end
    epsilon = (J1/Is)*mum/mue_si;
    Im = Iint/L;
    fe = fe_right/(L^2*Im);

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
    h = sqrt(mue*p);
    r = p/(1+e*cos(ta));
    n = sqrt(mue/a^3);
    u = aop + ta;

    % calculate force
    f_orbital = Iint*L*cross(uhat, [Bx, By, Bz]);
    a_orbital = f_orbital/m*1e-3'; % force to acceleraction in m to km
    a_r = a_orbital(1);
    a_h = a_orbital(2);
    a_theta = a_orbital(3);

    da = (2/(n*sqrt(1-e^2)))*(e*sin(ta)*a_r + (p/r)*a_theta);
    de = (sqrt(1-e^2)/(n*a))*(sin(ta)*a_r + ...
        (cos(ta)+(e+cos(ta))/(1+e*cos(ta)))*a_theta);
    di = (r*cos(u)/(n*a^2*sqrt(1-e^2)))*a_h;
    dRAAN = (r*sin(u)/(n*a^2*sqrt(1-e^2)*sin(i)))*a_h;
    daop = (sqrt(1-e^2)/(n*a*e))*(-cos(ta)*a_r + ...
        sin(ta)*(1+r/p)*a_theta) - (r*cot(i)*sin(u)/h)*a_r;
    dta = h/r^2 + (1/(e*h))*(p*cos(ta)*a_r - (p+r)*sin(ta)*a_theta);
    
    Tq = edt_torque_orbital(theta, phi, L, fe, Im, Bx, By, Bz);
    Tgg = gravity_grad_torque_orbital(r*1e3, theta, phi, Is, mue_si);
    T = Tq+Tgg;
    % x is towards zenith; 0
    % y is -h, rotation is in-plane or theta; Is
    % z is along track, rotation is out-of-plane or phi; Is
    ddtheta = T(2)/Is;
    ddphi = T(3)/Is;
%     ddtheta = 2*(1+dtheta)*dphi*tan(phi) - 1.5*sin(2*theta)...
%         -epsilon*(sin(i)*tan(phi)*(2*sin(ta)*cos(theta) - cos(ta)*sin(theta)) + cos(i));
%     ddphi = -sin(phi)*cos(phi)*((1+dtheta)^2 + 3*cos(theta)^2)...
%         +epsilon*sin(i)*(2*sin(ta)*sin(theta)+cos(ta)*cos(theta));
%     if t> 7.4
%         t
%     end
    dstate = [da; de; di; dRAAN; daop; dta; dtheta; dphi; ddtheta; ddphi];
    disp(dstate')
end
