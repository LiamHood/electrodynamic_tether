% function I = oml_sbet(h,p,L,At,PHI)
clear; close all; clc;
p = .0005*2*pi; %perimeter in meters
L = 5000;
At = .0005^2*pi;

me = 9.1093837e-31; % mass of electron
mi = 2.6561e-23;
Em = 165; % Induced electric field; dot(cross(v, B), u)
e = 1.60217663e-19; % electron charge
sigma = 3.5e7; % from paper
ht = 2*At/p; % characteristic transversal length
ninf = 2.0208e+11; % ionospheric plasma density 
mu = 1/172; % apparently for atomic oxygen
gamma = .15e-3;

delta = .5; %gamma*Em*L
Lstar = (me*Em)^(1/3)/(e*2^(7/3))*(3*pi*sigma*ht/ninf)^2/3;
lt = L/Lstar;
lt = 10;
% asymptotic expressions
i_B = mu*lt^(3/2)/2*(1+3/5*delta)*(1-mu^(2/3)*(3/2)*(1+delta)*(1+3/5*delta)^(-1/3)...
    -mu*(1/70)*lt^(3/2)*(81*delta^2+165*delta+70)/(5+3*delta));
xi_B = (2*i_B)^2/3*(1+4/15*i_B+53/360*i_B^2+3479/35640*i_B^3+3037/42768*i_B^4);

% curlye = h/Lstar;
% Isc = sigma*Em*At;
% I = Isc*ie(curlye);
% phie = PHI/(Em*Lstar);

% NEED TO SOLVE BVP

xi_low_mesh = linspace(0,xi_B,100);
solinit_low = bvpinit(xi_low_mesh,@guess,i_B);
sol_low = bvp4c(@bvp_ode_low, @bvp_bc_low, solinit_low);
xi_high_mesh = linspace(xi_B,lt,100);
solinit_high = bvpinit(xi_high_mesh,@guess,i_B);
sol_high = bvp4c(@bvp_ode_high, @bvp_bc_high, solinit_high);
figure
hold on
plot(sol_low.x,sol_low.y,'-o')
plot(sol_high.x,sol_high.y,'-o')

function dstate = bvp_ode_low(xi,state,i_B) % equation being solved
    ie = state(1);
    phie = state(2);
    
    dstate = zeros(2,1);
    dstate(1) = 3/4*sqrt(phie);
    dstate(2) = ie-1;
end

function res = bvp_bc_low(ya,yb,i_B) % boundary conditions    
res = [ya(1)
       yb(2)
       yb(1)-i_B];
end

function dstate = bvp_ode_high(xi,state,i_B) % equation being solved
    mu = 1/172;
    delta = .5;
    lt = 10;

    ie = state(1);
    phie = state(2);
    
    dstate = zeros(2,1);
    dstate(1) = -3/4*mu*sqrt(abs(phie))*(1+delta/lt*abs(phie));
    dstate(2) = ie-1;
end

function res = bvp_bc_high(ya,yb,i_B) % boundary conditions    
res = [yb(1)
       ya(2)
       ya(1)-i_B];
end

function y = guess(x) % guess at solution behavior
y = [1-x
     x^2];
end
