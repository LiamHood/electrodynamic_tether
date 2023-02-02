function Tq = edt_torque(states, tether_param, I, Bx, By, Bz)
    % assumes a uniform current flowing in the tether. In the case 
    % where the mass distribution is perfectly symmetrical, then the 
    % net torque about the center of mass is zero. This is an ideal 
    % scenario that is unlikely to be achieved in practice. This effect
    % will still be present in bare-wire tethers, but the relationship
    % is not as straightforward because of the nonuniform distribution 
    % of the electric current.
    roll = states(7);
    pitch = states(8);     

    % set states with friendly names
    L = tether_param(1)*1000;
    m1 = tether_param(2);
    m2 = tether_param(3);
    mt = tether_param(4);
    m = sum(tether_param(2:4));

    % nondimensional parameter that defines a measure of the average 
    % moment arm of the electromagnetic torque acting on the tether 
    % relative to the system center of mass
    PHI = (m1^2-m2^2+mt*(m1-m2))/m^2; 

    Qpitch = (I*L^2/2)*PHI*cos(roll)*...
        (sin(roll)*(Bx*cos(pitch)+By*sin(pitch))-Bz*cos(roll));
    Qroll = -(I*L^2/2)*PHI*(By*cos(pitch)-Bx*sin(pitch));
    Tq = [Qroll; Qpitch; 0]*1e-6; % convert to km
end
