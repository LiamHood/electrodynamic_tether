function [Bx, By, Bz] = MagField_OrbitalFrame(COES, mum)
    % Earthrs magnetic field vector in the Orbital Frame of SBET paper
    % SI units output

    % Define variables
    a = COES(1);
    i = COES(3);
    ta = COES(6);
    
    constant = mum/((a*1e3)^3);

    Bx = -2*constant*sin(i)*sin(ta);
    By = -constant*cos(i);
    Bz = constant*sin(i)*cos(ta);
end