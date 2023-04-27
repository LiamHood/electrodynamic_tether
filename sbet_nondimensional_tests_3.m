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

s0 = [0;1];

opts = optimoptions( 'fsolve' , 'Display' , 'off'  , 'FunctionTolerance' , 1e-8 , 'StepTolerance' , 1e-8 , 'OptimalityTolerance' , 1e-8 , 'Algorithm' , 'levenberg-marquardt' ) ; % 
[ s0s , F ] = fsolve( @TetherSolveFun , s0 , opts) ;
xi_span1 = [0, xi_B] ;
xi_span2 = [xi_B, lt] ;
opts = odeset( 'RelTol' , 1e-8 , 'AbsTol' , 1e-8 ) ;
[t1, s1] = ode45( @bvp_ode_low, xi_span1, s0s, opts) ;
sB = s1(end,:);
[t2, s2] = ode45( @bvp_ode_high, xi_span2, sB, opts) ;

figure
hold on
plot(t1, s1(:,1),'-o')
plot(t1, s1(:,2),'-o')
plot(t2, s2(:,1),'-o')
plot(t2, s2(:,2),'-o')
hold off
legend('i_e \xi<\xi_B', '\phi_e \xi<\xi_B', 'i_e \xi>\xi_B', '\phi_e \xi>\xi_B')

function dstate = bvp_ode_low(xi,state) % equation being solved
    ie = state(1);
    phie = state(2);
    
    dstate = zeros(2,1);
    dstate(1) = 3/4*sqrt(phie);
    dstate(2) = ie-1;
end

function dstate = bvp_ode_high(xi,state) % equation being solved
    mu = 1/172;
    delta = .5;
    lt = 10;

    ie = state(1);
    phie = state(2);
    
    dstate = zeros(2,1);
    dstate(1) = -3/4*mu*sqrt(abs(phie))*(1+delta/lt*abs(phie));
    dstate(2) = ie-1;
end

function F = TetherSolveFun(s0) % boundary conditions  
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

    xi_span1 = [0, xi_B] ;
    xi_span2 = [xi_B, lt] ;
    opts = odeset( 'RelTol' , 1e-8 , 'AbsTol' , 1e-8 ) ;
    [t1, s1] = ode45( @bvp_ode_low, xi_span1, s0, opts) ;
    sB = s1(end,:);
    [t2, s2] = ode45( @bvp_ode_high, xi_span2, sB, opts) ;
    sf = s2(end,:);
    % [ie, phie]

    F(1) = s0(1);
    F(2) = sB(1)-i_B;
    F(3) = sB(2);
    F(4) = sf(1);

end

function y = guess(x) % guess at solution behavior
    y = [1-x
         x^2];
end