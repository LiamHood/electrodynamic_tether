function Tgg = gravity_grad_torque_orbital(r, theta, phi, Is, mu)
    % Finds the gravity gradient torque effecting the tether
    % x is towards zenith
    % y is -h, rotation is in-plane or theta
    % z is along track, rotation is out-of-plane or phi
    inertia = [0, 0, 0; 0, Is, 0; 0, 0, Is];

    % Dynamics from "Spacecraft Dynamics and Control"
    u1 = [cos(phi)*cos(theta); -sin(phi); cos(phi)*sin(theta)];
    u2 = [-sin(theta); 0; cos(theta)];
    u3 = [-sin(phi)*cos(theta); -cos(phi); -sin(phi)*sin(theta)];
    tether_rotation = [u1, u2, u3];
    rb = tether_rotation*[r;0;0];

    % transform inertia into eci
    rbcross = [0, -rb(3), rb(2);
                rb(3), 0, -rb(1);
                -rb(2), rb(2), 0];
    Tgg = ((3*mu/r^5 * rbcross*inertia*rb));
end
    