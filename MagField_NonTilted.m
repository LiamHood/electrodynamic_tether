function [Bx, By, Bz] = MagField_NonTilted(COES)
    % Earth’s magnetic field vector in the Euler–Hill frame whose 
    % components are defined by (assuming a nontilted dipole)
    % Needs COES not states

    % Define variables
    a = COES(1);
    e = COES(2);
    i = COES(3);
    RAAN = COES(4);
    aop = COES(5);
    ta = COES(6);
    p = a*(1-e^2);
    
    % Calculate magnetic field based on https://doi.org/10.2514/1.12016
    Bo = 3.12e-5;
    Bo_over_R3 = Bo/((p/(1+e*cos(ta)))/6378)^3;

    Bx = -2*Bo_over_R3*sin(aop + ta)*sin(i);
    By = Bo_over_R3*cos(aop+ta)*sin(i);
    Bz = Bo_over_R3*cos(i);
end