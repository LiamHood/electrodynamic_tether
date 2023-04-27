% function I = oml_sbet(h,p,L,At,PHI)
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

parameters = [p, At, me, Em, e, sigma, ninf, mu, gamma];

% SOLVE BVP
I_B_0 = 1;
h0_0 = L/10;
s0 = [1; I_B_0; h0_0];
opts = optimoptions( 'fsolve' , 'Display' , 'off'  , 'FunctionTolerance' , 1e-8 ,...
    'StepTolerance' , 1e-8 , 'OptimalityTolerance' , 1e-8 , ...
    'Algorithm' , 'levenberg-marquardt') ; 
[ s0s , F ] = fsolve( @(s)TetherSolveFun(s,parameters,L) , s0 , opts) ;
h_span = [0, L];
parameters = [p, At, me, Em, e, sigma, ninf, mu, gamma, s0s(3)];
state0 = [0; s0s(1)];
opts = odeset( 'RelTol' , 1e-8 , 'AbsTol' , 1e-8 ) ;
[h, s] = ode45( @current_distribution_diffe, h_span, state0, opts, parameters) ;
fprintf('B is located at h0 = %f\n', s0s(3))
fprintf('Maximum current is I_B = %f\n', s0s(2))

J1 = 0;
for ii = 2:length(h)
    J1 = J1 + ((h_G-h(ii))*(s(ii,1))+(h_G-h(ii-1))*(s(ii-1,1)))/(h(ii)-h(ii-1));
end
fprintf('J1 = %f\n', J1)
figure
subplot(2,1,1)
plot(h, s(:,1),'-o')
xlabel('h [m]')
ylabel('I_e [A]')

subplot(2,1,2)
plot(h, s(:,2),'-o')
xlabel('h [m]')
ylabel('\Phi_e [V]')

