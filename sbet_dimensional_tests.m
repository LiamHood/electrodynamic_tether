% function I = oml_sbet(h,p,L,At,PHI)
clear; close all; clc;
L = 5000;
p = .0005*2*pi; %perimeter in meters
At = .0005^2*pi;
me = 9.1093837e-31; % mass of electron
Em = .165; % Induced electric field; dot(cross(v, B), u)
e = 1.60217663e-19; % electron charge
sigma = 3.5e7; % from paper
ninf = 2.0208e+11; % ionospheric plasma density 
mu = 1/172; % apparently for atomic oxygen
gamma = .15e-3;
parameters = [p, At, me, Em, e, sigma, ninf, mu, gamma];

% NEED TO SOLVE BVP
I_B_0 = 1;
h0_0 = L/10;
s0 = [0;1;I_B_0;h0_0];
opts = optimoptions( 'fsolve' , 'Display' , 'iter'  , 'FunctionTolerance' , 1e-8 ,...
    'StepTolerance' , 1e-8 , 'OptimalityTolerance' , 1e-8 , ...
    'Algorithm' , 'levenberg-marquardt','PlotFcn',@optimplotfirstorderopt) ; % 
[ s0s , F ] = fsolve( @TetherSolveFun , s0 , opts) ;
h_span1 = [0, s0s(4)] ;
h_span2 = [s0s(4), L] ;
opts = odeset( 'RelTol' , 1e-8 , 'AbsTol' , 1e-8 ) ;
[t1, s1] = ode45( @bvp_ode_low, h_span1, s0s(1:2), opts, parameters) ;
sB = s1(end,:);
[t2, s2] = ode45( @bvp_ode_high, h_span2, sB, opts, parameters) ;
sf = s2(end,:);
[tL, sL] = ode45( @current_distribution_diffe, [0, L], s0s(1:2), opts, [parameters, s0s(4)]) ;
fprintf('B location initial guess is h0 = %f\n', h0_0)
fprintf('Maximum current initial guess is I_B = %f\n', I_B_0)
fprintf('B is located at h0 = %f\n', s0s(4))
fprintf('Maximum current is I_B = %f\n', s0s(3))
fprintf('PHI_e at A is %f\n', s0s(2))
fprintf('PHI_e at B is %f\n', sB(2))
fprintf('PHI_e at C is %f\n', sf(2))
fprintf('I_e at A is %f\n', s0s(1))
fprintf('I_e at B is %f\n', sB(1))
fprintf('I_e at C is %f\n', sf(1))
figure
subplot(2,1,1)
hold on
plot(t1, s1(:,1),'-o')
plot(t2, s2(:,1),'-o')
plot(tL, sL(:,1),'k')
hold off
legend('I_e h<h_0', 'I_e h>h_0', 'I_e')
xlabel('h')
ylabel('I_e')

subplot(2,1,2)
hold on
plot(t1, s1(:,2),'-o')
plot(t2, s2(:,2),'-o')
plot(tL, sL(:,2),'k')
hold off
legend('\Phi_e h<h_0', '\Phi_e h>h_0', '\Phi_e')
xlabel('h')
ylabel('\Phi_e')

function dstate = current_distribution_diffe(h,state,parameters) % equation being solved

    Ie = state(1);
    PHIe = state(2);

    p = parameters(1); %perimeter in meters
    At = parameters(2);
    me = parameters(3); % mass of electron
    Em = parameters(4); % Induced electric field; dot(cross(v, B), u)
    e = parameters(5); % electron charge
    sigma = parameters(6); % from paper
    ninf = parameters(7); % ionospheric plasma density 
    mu = parameters(8); % apparently for atomic oxygen
    gamma = parameters(9);
    h0 = parameters(10);

    dstate = zeros(2,1);
    if real(h) < h0
        dstate(1) = e*ninf*p/pi*sqrt((2*e/me)*PHIe);
    else
        dstate(1) = -e*ninf*p/pi*sqrt((2*e/me)*abs(PHIe))*mu*(1+gamma*abs(PHIe));
    end
    dstate(2) = Ie/(sigma*At)-Em;
end

function dstate = bvp_ode_low(h,state,parameters) % equation being solved
    
    Ie = state(1);
    PHIe = state(2);
    
    p = parameters(1); %perimeter in meters
    At = parameters(2);
    me = parameters(3); % mass of electron
    Em = parameters(4); % Induced electric field; dot(cross(v, B), u)
    e = parameters(5); % electron charge
    sigma = parameters(6); % from paper
    ninf = parameters(7); % ionospheric plasma density 
    mu = parameters(8); % apparently for atomic oxygen
    gamma = parameters(9);

    dstate = zeros(2,1);
    dstate(1) = e*ninf*p/pi*sqrt((2*e/me)*PHIe);
    dstate(2) = Ie/(sigma*At)-Em;
end

function dstate = bvp_ode_high(h,state,parameters) % equation being solved

    Ie = state(1);
    PHIe = state(2);

    p = parameters(1); %perimeter in meters
    At = parameters(2);
    me = parameters(3); % mass of electron
    Em = parameters(4); % Induced electric field; dot(cross(v, B), u)
    e = parameters(5); % electron charge
    sigma = parameters(6); % from paper
    ninf = parameters(7); % ionospheric plasma density 
    mu = parameters(8); % apparently for atomic oxygen
    gamma = parameters(9);

    dstate = zeros(2,1);
    dstate(1) = -e*ninf*p/pi*sqrt((2*e/me)*abs(PHIe))*mu*(1+gamma*abs(PHIe));
    dstate(2) = Ie/(sigma*At)-Em;
end

function F = TetherSolveFun(s0) % boundary conditions  
    p = .0005*2*pi; %perimeter in meters
    At = .0005^2*pi;
    me = 9.1093837e-31; % mass of electron
    Em = .165; % Induced electric field; dot(cross(v, B), u)
    e = 1.60217663e-19; % electron charge
    sigma = 3.5e7; % from paper
    ninf = 2.0208e+11; % ionospheric plasma density 
    mu = 1/172; % apparently for atomic oxygen
    gamma = .15e-3;
    parameters = [p, At, me, Em, e, sigma, ninf, mu, gamma];
    L = 5000;

    I_B = s0(3);
    h0 = s0(4);
    h_span1 = [0, h0] ;
    h_span2 = [h0, L] ;
    opts = odeset( 'RelTol' , 1e-8 , 'AbsTol' , 1e-8 ) ;
    [t1, s1] = ode45( @bvp_ode_low, h_span1, s0(1:2), opts, parameters) ;
    sB = s1(end,:);
    [t2, s2] = ode45( @bvp_ode_high, h_span2, sB, opts, parameters) ;
    sf = s2(end,:);
    % [ie, phie]

    F(1) = s0(1)*10;
    F(2) = sB(1)-I_B;
    F(3) = sB(2);
    F(4) = sf(1)*10;
    F(5) = imag(s0(1));
    F(6) = imag(sB(1));
    F(7) = imag(sf(1));
    F(8) = imag(s0(2));
    F(9) = imag(sB(2));
    F(10) = imag(sf(2));
    F(11) = imag(s0(3));
    F(12) = imag(s0(4));
end
