clear; close all; clc;

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

% parameters = [p, At, me, Em, e, sigma, ninf, mu, gamma];

