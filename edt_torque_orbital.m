function Tq = edt_torque_orbital(theta, phi, L, fe, Im, Bx, By, Bz)
    % x is towards zenith
    % y is -h, rotation is in-plane or theta
    % z is along track, rotation is out-of-plane or phi

    Qtheta = fe*Im*L^2*cos(phi)*(sin(phi)*(Bx*cos(theta)+Bz*sin(theta))+By*cos(phi));
    Qphi = -fe*Im*L^2*(Bx*sin(theta)-Bz*cos(theta));
    Tq = [0; Qtheta; Qphi];
end
