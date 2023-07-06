function [Bx, By, Bz] = MagField_OrbitalFrame(COES)
    % Earthrs magnetic field vector in the Orbital Frame of SBET paper
    % SI units output

    % Define variables
    a = COES(1);
    i = COES(3);
    ta = COES(6);
    
%     Bo = 3.12e-5;
%     Bo_over_R3 = Bo/((p/(1+e*cos(ta)))/6378)^3;
    mum = 8.22e22; % https://www.britannica.com/science/geomagnetic-field/Dipolar-field
    constant = mum/((a*1e3)^3);

    Bx = -2*constant*sin(i)*sin(ta);
    By = -constant*cos(i);
    Bz = constant*sin(i)*cos(ta);
end