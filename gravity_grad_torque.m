function Tgg = gravity_grad_torque(states, tether_param, mu)
    % Finds the gravity gradient torque effecting the tether

    % set states with friendly names
    a = states(1);
    e = states(2);
    i = states(3);
    RAAN = states(4);
    aop = states(5);
    ta = states(6);
    roll = states(7);
    pitch = states(8);
    yaw = states(9);

    % Make inertia matrix
    Ix = tether_param(5);
    Iy = tether_param(6);
    Iz = tether_param(7);
    inertia = [Ix, 0, 0; 0, Iy, 0; 0, 0, Iz];
    
    % Calculate the position and velocity vector
    [ rvec , vvec ] = coes2state( [sqrt(mu*a*(1-e^2)), i, e, RAAN, aop, ta] , mu );
    r = norm(rvec);

    % find Roll-Pitch-Yaw axes
%     yaw_hat = -rvec/r;
%     pitch_hat = cross(rvec,vvec)/norm(cross(rvec,vvec));
%     roll_hat = cross(yaw_hat, pitch_hat)/norm(cross(yaw_hat, pitch_hat));
%     C_rpy2eci = [roll_hat, pitch_hat, yaw_hat];
%     C_eci2rpy = C_rpy2eci';

    % Dynamics from "Spacecraft Dynamics and Control"
    C_rpy2body = [cos(pitch)*cos(yaw), cos(pitch)*sin(yaw), -sin(pitch);
            sin(roll)*sin(pitch)*cos(yaw)-cos(roll)*sin(yaw), ...
            sin(roll)*sin(pitch)*sin(yaw)+cos(roll)*cos(yaw), ...
            sin(roll)*cos(pitch);
            cos(roll)*sin(pitch)*cos(yaw)+sin(roll)*sin(yaw), ...
            cos(roll)*sin(pitch)*sin(yaw)-sin(roll)*cos(yaw), ...
            cos(roll)*cos(pitch)];
    C_body2rpy = C_rpy2body';
    rb = C_rpy2body*[0;0;r];

    % transform inertia into eci
    rbcross = [0, -rb(3), rb(2);
                rb(3), 0, -rb(1);
                -rb(2), rb(2), 0];
    Tgg = ((3*mu/r^5 * rbcross*inertia*rb));
end
    