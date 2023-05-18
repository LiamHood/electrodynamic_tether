clear; close all; clc;
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
i = deg2rad(60); % [degree], inclination
RAAN = deg2rad(0); % [degree], right ascension of ascending node
aop = deg2rad(0); % [degree], argument of perigee
ta = deg2rad(0); % [degree], true anomaly
coes_state0 = [a; e; i; RAAN; aop; ta];
[r0, v0] = coes2state([sqrt(mue*a*(1-e^2)),i,e,RAAN,aop,ta],mue);

tol = 1e-8;
tspan = [0, 6*60*60];
[ t , r , v ] = TwoBody( tspan , r0 , v0 , mue , tol );
for ii = 1:length(t)
    ihato_m = r(:,ii)/norm(r(:,ii));
    jhato_m = cross(r(:,ii),v(:,ii))/norm(cross(r(:,ii),v(:,ii)));
    khato_m = cross(ihato_m,jhato_m);
    orbital2moving = [ihato_m, jhato_m, khato_m];
    uhat_m = orbital2moving*[1;0;0];
    v_o = orbital2moving'*v(:,ii);
    COES = state2coes([r(:,ii);v(:,ii)], mue);
    [Bx, By, Bz] = MagField_NonTilted(COES);
    Em_first = cross(v_o*1e3, [Bx; By; Bz]);
    Em(ii) = dot(cross(v_o*1e3, [Bx; By; Bz]), [1;0;0]);
    rmag(ii) = norm(r(:,ii));
end


figure
plot(t, Em)


figure
plot(rmag, Em)