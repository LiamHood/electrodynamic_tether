function [fr, ftheta, fh] = edt_forces(states, tether_param, I,  Bx, By, Bz)
    % assumes a uniform current flowing in the tether
    
    % set states with friendly names
    theta = states(8);
    phi = states(7);

    L = tether_param(1)*1000; % change to meters
    m = sum(tether_param(2:4));

    fr = (I*L/m*(Bz*sin(theta)*cos(phi)-By*sin(phi)))*1e-3; 
    ftheta = (I*L/m*(Bx*sin(phi)-Bz*cos(theta)*cos(phi)))*1e-3;
    fh = (I*L/m*(By*cos(theta)*cos(phi)-Bx*sin(theta)*cos(phi)))*1e-3;
end
