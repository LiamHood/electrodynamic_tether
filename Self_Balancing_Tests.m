clear; close all; clc;
% https://arc-aiaa-org.proxy.wichita.edu/doi/pdf/10.2514/6.2004-5309
% Self Balancing Tether
% m2 is upper mass
% m1 is lower mass
% mt is tether mass
% (m1, m2, mt) => (m, phi, LAMBDAt)

m1 = 100;
m2 = 10;
mt = 10;
L = 5000;
[m, phi, LAMBDAt, hG, Is] = params_2_sbet_params(m1, m2, mt, L);
%rotate about y then about x'
% x is radial, y is momentum, z is 'velocity'
a = 7378; % [km], semi-major axis
e = .002; % eccentricity
i = deg2rad(30); % [degree], inclination
RAAN = deg2rad(40); % [degree], right ascension of ascending node
aop = deg2rad(50); % [degree], argument of perigee
ta = deg2rad(0); % [degree], true anomaly
mue = 398600; % mu of earth
mum = 3.12e-5; % magnetic constant of earths mag field
[r, v] = coes2state([sqrt(mue*a*(1-e^2)),i,e,RAAN,aop,ta],mue);

ihati = r/norm(r);
jhati = cross(r,v)/norm(cross(r,v));
khati = cross(ihati,jhati);
orbital2inertial = [ihati, jhati, khati];
ihato = [1;0;0] ; % unit vector of orbit frame
jhato = [0;1;0];
khato = [0;0;1];

in_plane = deg2rad(10);
out_plane = deg2rad(5);
% [cos(out_plane), 0, sin(out_plane); 0, 1, 0; -sin(out_plane), 0, cos(out_plane)]
% [1, 0, 0; 0, cos(out_plane), -sin(out_plane); 0, sin(out_plane), cos(out_plane)]
orbital2body = [cos(out_plane), 0, sin(out_plane); 0, 1, 0; -sin(out_plane), 0, cos(out_plane)]*...
    [cos(in_plane), -sin(in_plane), 0; sin(in_plane), cos(in_plane), 0; 0, 0, 1];
% orbital2body = [cos(in_plane), -sin(in_plane), 0; sin(in_plane), cos(in_plane), 0; 0, 0, 1]*...
%     [cos(out_plane), 0, sin(out_plane); 0, 1, 0; -sin(out_plane), 0, cos(out_plane)];
% uhat0 = orbital2body'*ihato % unit vector of tether
uhat = cos(out_plane)*cos(in_plane)*ihato - sin(out_plane)*jhato + cos(out_plane)*sin(in_plane)*khato;
wr = sqrt(mu/a^3);

syms h
J1 = integral((hG-h)*Ie(h),0,L);
epsilon = (J1/Is)*(mum/mue);

tol = 1e-10;
opts = odeset('RelTol', tol, 'AbsTol', tol ) ;
[ t , states ] = ode45(@sbet_prop, tspan , ...
    states0, opts, epsilon, mue) ;

function dstate = sbet_prop(t, states, epsilon, mu)
    a = states(1);
    e = states(2);
    i = states(3);
    RAAN = states(4);
    aop = states(5);
    ta = states(6);

    in_plane = states(7);
    out_plane = states(8);
    dip = states(9);
    dop = states(10);

    fr = 0;
    fs = 0;
    fw = 0;
    da = (2/(n*sqrt(1-e^2)))*(e*sin(ta)*fr + (p/r)*fs);
    de = (sqrt(1-e^2)/(n*a))*(sin(ta)*fr + ...
        (cos(ta)+(e+cos(ta))/(1+e*cos(ta)))*fs);
    di = (r*cos(u)/(n*a^2*sqrt(1-e^2)))*fw;
    dRAAN = (r*sin(u)/(n*a^2*sqrt(1-e^2)*sin(i)))*fw;
    daop = (sqrt(1-e^2)/(n*a*e))*(-cos(ta)*fr + ...
        sin(ta)*(1+r/p)*fs) - (r*cot(i)*sin(u)/h)*fw;
    dta = h/r^2 + (1/(e*h))*(p*cos(ta))*fr - (p+r)*sin(ta)*fs;

    ddip = 2*(1+dip)*dop*tan(out_plane) - 3/2*sin(2*in_plane) - epsilon*(sin(i)*tan(out_plane)*(2*sin(ta)*cos(in_plane)-cos(ta)*cos(in_plane))+cos(i));
    ddop = -sin(out_plane)*cos(out_plane)*((1+dip)^2+3*cos(in_plane)^2) + epsilon*sin(i)*(2*sin(ta)*sin(in_plane)+cos(ta)*cos(in_plane));

    dstate = [da, de, di, dRAAN, daop, dta, dip, dop, ddip, ddop];
end

function I = oml_sbet(h,e,ninf,p,me,Em,sigma,ht,L,At,PHI)
    ht = 2*At/p;
    Lstar = (me*Em)^(1/3)/(e*2^(7/3))*(3*pi*sigma*ht/ninf)^2/3;
    curlye = h/Lstar;
    lt = L/Lstar;
    Isc = sigma*Em*At;
    phie = PHI/(Em*Lstar);
    I = Isc*ie(curlye);

end

% % Gravitational Torque
% Mg = 3*mue/a^3*Is*cross(uhat,ihato)*dot(uhat,ihato)
% 
% syms h
% Bx = -2*mum/a^3*sin(i)*sin(ta);
% By = -mum/a^3*cos(i);
% Bz = mum/a^3*sin(i)*cos(ta);
% J1 = integral((hG-h)*Ie(h),0,L);
% Me = cross(uhat, cross(uhat, B))*J1
