clear; close all; clc;
% https://arc-aiaa-org.proxy.wichita.edu/doi/pdf/10.2514/6.2004-5309

%% Define tether parameters
L = 5000; % tether length
tether_radius = .0005;
per = tether_radius*2*pi; %perimeter in meters
At = tether_radius^2*pi;

%% Define constants
mue = 398600; % mu of earth
mue_si = mue*1e9; % mu of earth in meters
mum = 8e15; % magnetic constant of earths mag field
me = 9.1093837e-31; % mass of electron
mi = 2.6561e-26; % mass of atomic oxygen
e = 1.60217663e-19; % electron charge
mu = sqrt(me/mi); % Ratio (me/mi)^1/2
sigma = 3.5e7; % conductivity for aluminum
gamma = .15e-3; % secondary emission yield

%% Define variables
nn = 1e2;
mm = 1e2;
Em = linspace(0,1e6,nn);
ninf = linspace(2e9,2e12,mm);

for ii = 1:nn
    for jj = 1:mm
        percent_complete = 100*((ii-1)*mm+jj)/(nn*mm)
        parameters = [per, At, me, Em(ii), e, sigma, ninf(jj), mu, gamma];
        % SOLVE BVP
        I_B_0 = 1;
        h0_0 = L/10;
        s0 = [1;I_B_0;h0_0];
        opts = optimoptions( 'fsolve' , 'Display' , 'none', 'FunctionTolerance' , 1e-8 ,...
            'StepTolerance' , 1e-8 , 'OptimalityTolerance' , 1e-8 , ...
            'Algorithm' , 'levenberg-marquardt') ; 
        [ s0s , F, exitflag_iter] = fsolve( @(s)TetherSolveFun(s,parameters,L) , s0 , opts) ;
        PHIe_0(ii,jj) = s0s(1); 
        I_B(ii,jj) = s0s(2);
        h_0(ii,jj) = s0s(3);
        exitflag(ii,jj) = exitflag_iter;
%         parameters = [per, At, me, Em(ii), e, sigma, ninf(jj), mu, gamma, h_0];
%         opts = odeset( 'RelTol' , 1e-8 , 'AbsTol' , 1e-8 ) ;
%         [h, s] = ode45( @current_distribution_diffe, [0, L], [0, PHIe_0], opts, parameters) ;
%         Iint = 0;
%         for ii = 2:length(h)
%             Iint = Iint + (s(ii,1)+s(ii-1,1))/2*(h(ii)-h(ii-1));
%         end
    end
end
save("starting_values_1e2_1e2_1e6_5km_5mm.mat","Em", "ninf", "PHIe_0", "I_B", "h_0", "exitflag")
figure
surf(ninf, Em, exitflag)