clear; close all; clc;
% https://arc-aiaa-org.proxy.wichita.edu/doi/pdf/10.2514/6.2004-5309

%% Define tether
tether.m1 = 321.6; % lower mass [kg]
tether.m2 = 563.5; % upper mass [kg]
tether.L = 20000; % tether length [km]
tether.density = 2700; % [kg/m^3]
tether.thickness = .18e-3; % [m]
tether.perimeter = 24e-3; % [m]
tether.At = (tether.perimeter-tether.thickness*2)/2*tether.thickness; % [m^2]
tether.mt = tether.density*tether.L*tether.At; % tether mass [kg]
[tether.m, tether.phi, tether.LAMBDAt, tether.h_G, tether.Is] = params_2_sbet_params(tether.m1, tether.m2, tether.mt, tether.L); % [kg, radians, non-dimensional, m, kg*m^2]
tether.sigma = 3.5e7; % conductivity for aluminum [1/(Ohm*m)]

%% define constants
constants.mue = 398600; % mu of earth [km^3/s^2]
constants.mue_si = constants.mue*1e9; % mu of earth in meters [m^3/s^2]
A_m2_mum = 8.22e22; % magnetic constant of earths mag field [A*m^2] https://www.britannica.com/science/geomagnetic-field/Dipolar-field
constants.mum = A_m2_mum/795774.71545947673925; % magnetic constant of earths mag field [T*M^3] https://maurermagnetic.com/en/demagnetizing/technology/convert-magnetic-units/
constants.me = 9.1093837e-31; % mass of electron [kg]
constants.mi = 2.6561e-26; % mass of atomic oxygen [kg]
constants.e = 1.60217663e-19; % electron charge [C]
constants.mu = sqrt(constants.me/constants.mi); % Ratio 
constants.gamma = .15e-3; % secondary emission yield [1/V]

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
Em = linspace(1,2,nn);
ninf = linspace(2e9,2e12,mm);

for ii = 1:nn
    for jj = 1:mm
        parameters.Em = Em(ii);
        parameters.ninf = ninf(jj);
        percent_complete = 100*((ii-1)*mm+jj)/(nn*mm)
        % SOLVE BVP
        I_B_0 = 1;
        h0_0 = tether.L/10;
        s0 = [1;I_B_0;h0_0];
        opts = optimoptions( 'fsolve' , 'Display' , 'none'  , 'FunctionTolerance' , 1e-8 ,...
            'StepTolerance' , 1e-8 , 'OptimalityTolerance' , 1e-8 , ...
            'Algorithm' , 'levenberg-marquardt') ; 
        [ s0s , F , exitflag_iter] = fsolve(@(s)TetherSolveFun2(s, tether, constants, parameters), s0 , opts) ;
        PHIe_0(ii,jj) = s0s(1); 
        h_0(ii,jj) = s0s(2);
        exitflag(ii,jj) = exitflag_iter;

    end
end
save("Em_1_2_1e2_ninf_2e9_2e12_1e2.mat","Em", "ninf", "PHIe_0", "h_0", "exitflag")
figure
surf(ninf, Em, exitflag)