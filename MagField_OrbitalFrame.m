function [Bx, By, Bz] = MagField_NonTilted(COES)
    % Earthrs magnetic field vector in the Orbital Frame of SBET paper
    % SI units output

    % Define variables
    a = COES(1);
    e = COES(2);
    i = COES(3);
    RAAN = COES(4);
    aop = COES(5);
    ta = COES(6);
    p = a*(1-e^2);
    
%     Bo = 3.12e-5;
%     Bo_over_R3 = Bo/((p/(1+e*cos(ta)))/6378)^3;
    mum = 8e15; % magnetic constant of earths mag field
    constant = mum/((a*1e3)^3);

    Bx = -2*constant*sin(i)*sin(ta);
    By = -constant*cos(i);
    Bz = constant*sin(i)*cos(ta);
end