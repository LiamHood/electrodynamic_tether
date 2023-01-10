function Tgg = gravity_grad_torque(states, tether_param, mu)
    % Finds the gravity gradient torque effecting the tether

    % set states with friendly names
    a = states(1);
    e = states(2);
    i = states(3);
    RAAN = states(4);
    aop = states(5);
    ta = states(6);
    phi = states(7);
    theta = states(8);
    psi = states(9);

    % Make inertia matrix
    Ix = tether_param(5);
    Iy = tether_param(6);
    Iz = tether_param(7);
    inertia = [Ix, 0, 0; 0, Iy, 0; 0, 0, Iz];
    
    % Calculate the position and velocity vector
    [ rvec , vvec ] = coes2state( [sqrt(mu*a*(1-e^2)), i, e, RAAN, aop, ta] , mu );
    r = norm(rvec);

    % Dynamics from "Spacecraft Dynamics and Control"
    Cbo = [cos(theta)*cos(psi), cos(theta)*sin(psi), -sin(theta);
            sin(phi)*sin(theta)*cos(psi)-cos(theta)*sin(psi), ...
            sin(phi)*sin(theta)*sin(psi)+cos(phi)*sin(psi), ...
            sin(phi)*cos(theta);
            cos(phi)*sin(theta)*cos(psi)+sin(phi)*sin(psi), ...
            cos(phi)*sin(theta)*sin(psi)-sin(phi)*cos(psi), ...
            cos(phi)*cos(theta)];
    rb = Cbo*[0;0;-r];

    % transform inertia into eci
    rbcross = [0, -rb(3), rb(2);
                rb(3), 0, -rb(1);
                -rb(2), rb(2), 0];
    Tgg = 3*mu/r^5 * rbcross*inertia*rb;
end
    