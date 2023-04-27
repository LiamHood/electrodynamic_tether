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
syms h
J1inside = @(h) (h_G-h)*Ie(h);
J1 = integral(J1inside,0,L);
epsilon = (J1/Is)*mum/mue_si;


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
% ihato = [1;0;0] ; % unit vector of orbit frame
% jhato = [0;1;0];
% khato = [0;0;1];

in_plane = deg2rad(10); % theta
out_plane = deg2rad(5); % phi
uhato = [cos(out_plane)*cos(in_plane); -sin(out_plane); cos(out_plane)*sin(in_plane)];
% [MG, ME] = sbet_torque(mue, a, i, ta, Is, orbital2moving*uhato, h_G, L);
tol = 1e-10;
tspan = [0, 200];
tether_state0 = [in_plane; out_plane; 0; 0];
opts = odeset('RelTol', tol, 'AbsTol', tol ) ;
[ t , states ] = ode45(@sbet_libration, tspan , [coes_state0; tether_state0], opts, mue, epsilon) ;

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

% 

% [cos(out_plane), 0, sin(out_plane); 0, 1, 0; -sin(out_plane), 0, cos(out_plane)]
% [1, 0, 0; 0, cos(out_plane), -sin(out_plane); 0, sin(out_plane), cos(out_plane)]
% orbital2body = [cos(out_plane), 0, sin(out_plane); 0, 1, 0; -sin(out_plane), 0, cos(out_plane)]*...
%     [cos(in_plane), -sin(in_plane), 0; sin(in_plane), cos(in_plane), 0; 0, 0, 1];
% orbital2body = [cos(in_plane), -sin(in_plane), 0; sin(in_plane), cos(in_plane), 0; 0, 0, 1]*...
%     [cos(out_plane), 0, sin(out_plane); 0, 1, 0; -sin(out_plane), 0, cos(out_plane)];
% uhat0 = orbital2body'*ihato % unit vector of tether
% uhat = cos(out_plane)*cos(in_plane)*ihato - sin(out_plane)*jhato + cos(out_plane)*sin(in_plane)*khato;
% wr = sqrt(mue/a^3);
% 
% syms h
% J1 = integral((hG-h)*Ie(h),0,L);
% epsilon = (J1/Is)*(mum/mue);
% 
% tol = 1e-10;
% opts = odeset('RelTol', tol, 'AbsTol', tol ) ;
% [ t , states ] = ode45(@sbet_prop, tspan , ...
%     states0, opts, epsilon, mue) ;
% 


function I = oml_sbet(h,p,L,At,PHI)
    me = 9.1093837e-31; % mass of electron
    mi = 2.6561e-26;
    Em = 165; % Induced electric field; dot(cross(v, B), u)
    e = 1.60217663e-19; % electron charge
    sigma = 3.5e7; % from paper
    ht = 2*At/p; % characteristic transversal length
    ninf = 2.0208e+11; % ionospheric plasma density 
    mu = 1/172; % apparently for atomic oxygen
    gamma = 1;
    
    delta = .5; %gamma*Em*L
    Lstar = (me*Em)^(1/3)/(e*2^(7/3))*(3*pi*sigma*ht/ninf)^2/3;
    lt = L/Lstar;
    % asymptotic expressions
    i_B = mu*lt^(3/2)/2*(1+3/5*delta)*(1-mu^(2/3)*(3/2)*(1+delta)*(1+3/5*delta)^(-1/3)...
        -mu*(1/70)*lt^(3/2)*(81*delta^2+165*delta+70)/(5+3*delta));
    curlye_B = (2*i_B)^2/3*(1+4/15*i_B+53/360*i_B^2+3479/35640*i_B^3+3037/42768*i_B^4);
    
    curlye = h/Lstar;
    Isc = sigma*Em*At;
    I = Isc*ie(curlye);
    phie = PHI/(Em*Lstar);
    
    % NEED TO SOLVE BVP

end
