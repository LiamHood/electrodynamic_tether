function [Bx, By, Bz] = MagField_igrf(states, t, mu)
    % Earth’s magnetic field vector in the Euler–Hill frame whose 
    % components are defined by (assuming a nontilted dipole)

    % Define variables
    a = states(1);
    e = states(2);
    i = states(3);
    RAAN = states(4);
    aop = states(5);
    ta = states(6);
    p = a*(1-e^2);
    
    % calculate position a few ways
    [ rvec , vvec ] = coes2state( [sqrt(mu*a*(1-e^2)), i, e, RAAN, aop, ta] , mu );
    tsec = t;
    tmin = 0;
    thour = 12;
    tday = 1;
    while tsec >= 60
        tsec = tsec - 60;
        tmin = tmin + 1;
    end
    while tmin >= 60
        tmin = tmin - 60;
        thour = thour + 1;
    end
    while thour >= 24
        thour = thour - 24;
        tday = tday + 1;
    end
    time = [2020 1 tday thour tmin tsec];
    [ rvec_ecef , vvec_ecef ] = eci2ecef(time, rvec, vvec);
    LLA = ecef2lla(rvec_ecef'*1000);
%         LLA = eci2lla(rvec'.*1000,time);
    if LLA(3) < 0
        LLA(3) = (norm(rvec)-6378)*1000;
    end
    [B_NED,~,~,~,~] = igrfmagm(LLA(3),LLA(1),LLA(2),decyear(time),13);
    [B_ECEFX, B_ECEFY, B_ECEFZ] = ned2ecefv(B_NED(1),B_NED(2),B_NED(3),LLA(1),LLA(2));
    [r_eci, B_eci] = ecef2eci(time, rvec_ecef, [B_ECEFX, B_ECEFY, B_ECEFZ]);
    y_eh_ineci = r_eci/norm(r_eci);
    z_eh_ineci = cross(rvec,vvec)/norm(cross(rvec,vvec));
    x_eh_ineci = cross(y_eh_ineci,z_eh_ineci)/norm(cross(y_eh_ineci,z_eh_ineci));
    eci2eh = [x_eh_ineci, y_eh_ineci, z_eh_ineci];
    B_eh = eci2eh*B_eci;

    r_rsw = rvec/norm(rvec);
    w_rsw = cross(rvec,vvec)/norm(cross(rvec,vvec));
    s_rsw = cross(w_rsw, r_rsw)/norm(cross(w_rsw, r_rsw));
    rsw2eci = [r_rsw, s_rsw, w_rsw];

    %NEED IN LVLH
    B_eh = rsw2eci'*B_eh*1e-9;
    Bx = B_eh(1);
    By = B_eh(2);
    Bz = B_eh(3);

end
