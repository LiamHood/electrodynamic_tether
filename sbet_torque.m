function [MG, ME] = sbet_torque(mu, a, i, ta, Is, uhat_m, h_G, L)
    % Find torques for describing motion relative to the Gx1y1z1 moving
    % frame, with origin at center of mass of tether and axes parallel to
    % the inertial frame, ECI
    % mu is gravitational constant of earth
    % a is semi-major axis
    % uhat_m is unit vector of tether in the moving frame
    % ihat1_m is the unit vector defining x of moving frame in moving frame
    ihat1_m = [1; 0; 0];
    MG = 3*mu/(a^3) * Is * cross(uhat_m, ihat1_m) * dot(uhat_m, ihat1_m);
    syms h
    J1inside = @(h) (h_G-h)*Ie(h);
    J1 = integral(J1inside,0,L);
    B = sbet_magfield(a, i, ta);
    ME = cross(uhat_m, cross(uhat_m, B))*J1;
end