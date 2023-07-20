function Tgg = gravity_grad_torque_orbital(r, theta, phi, Is, mu)
    % Finds the gravity gradient torque effecting the tether
    % x is towards zenith
    % y is -h, rotation is in-plane or theta
    % z is along track, rotation is out-of-plane or phi
    inertia = [0, 0, 0; 0, Is, 0; 0, 0, Is];

    C_attitude = tether_attitude(theta, phi);
    rb = C_attitude*[r;0;0];

    % transform inertia into orbital
    rbcross = [0, -rb(3), rb(2);
                rb(3), 0, -rb(1);
                -rb(2), rb(2), 0];
    Tgg = ((3*mu/r^5 * rbcross*inertia*rb));
end
    