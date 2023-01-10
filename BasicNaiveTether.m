function [ t , states] = BasicNaiveTether( tspan , sc_state0, tether_state0, tether_param, mu , tol )
%     Based on https://arc.aiaa.org/doi/10.2514/1.12016
% Uses Classical Orbital Elements a (semimajor axis), e (eccentricity), 
% i (inclination), RAAN (right ascension of ascending node), aop (argument
% of periapsis), ta (true anomaly)
%
% Considers electrodynamic force,
% will at some point maybe: atmospheric drag, spherical harmonics,
% third body of sun and moon, srp
% tether is treated as dumbbell model, rigid with lumped masses at the end

% tether_param = L, m1, m2, mt

    opts = odeset('RelTol', tol, 'AbsTol', tol ) ;
    [ t , states ] = ode45(@gauss_variations, tspan , ...
        [ sc_state0 ; tether_state0], opts, tether_param, mu) ;


    function dstate = gauss_variations(t, states, tether_param, mu)
        a = states(1);
        e = states(2);
        i = states(3);
        RAAN = states(4);
        aop = states(5);
        ta = states(6);

        theta = states(7);
        phi = states(8);
        dtheta = states(9);
        dphi = states(10);

        I = states(11);

        L = tether_param(1);
        m1 = tether_param(2);
        m2 = tether_param(3);
        mt = tether_param(4);
        m = sum(tether_param(2:4));

        [fr, fs, fw] = edt_forces(states, tether_param);
        [Qtheta, Qphi] = edt_torque(states, tether_param);
        [GGtheta, GGphi] = gravity_grad_torque(states, tether_param, mu);

        p = a*(1-e^2);
        h = sqrt(mu*p);
        r = norm(states(1:3));
        n = sqrt(mu/a^3);
        u = aop + ta;
        mstar = (m1+mt/2)*(m2+mt/2)/m-mt/6;

        % orbital elements using Vallado p.636
        da = (2/(n*sqrt(1-e^2)))*(e*sin(ta)*fr + (p/r)*fs);
        de = (sqrt(1-e^2)/(n*a))*(sin(ta)*fr + ...
            (cos(ta)+(e+cos(ta))/(1+e*cos(ta)))*fs);
        di = (r*cos(u)/(n*a^2*sqrt(1-e^2)))*fw;
        dRAAN = (r*sin(u)/(n*a^2*sqrt(1-e^2)*sin(i)))*fw;
        daop = (sqrt(1-e^2)/(n*a*e))*(-cos(ta)*fr + ...
            sin(ta)*(1+r/p)*fs) - (r*cot(i)*sin(u)/h)*fw;
        dta = h/r^2 + (1/(e*h))*(p*cos(ta))*fr - (p+r)*sin(ta)*fs;
        
        % second derivative of true anomaly is found by assuming that most
        % of the change in true anomaly is due to the h/r^2 term and then
        % uses 9-19 from Vallado for dh/dt, (dh/dt = r*Fs)
%         ddta = fs/r;
        ddta = 0;

        % libration dynamics from Paul Williams paper
        ddtheta = -ddta + 2*(dtheta + dta)*dphi*tan(phi) ...
            - 3*(dta^2/(1+e*cos(ta)))*sin(theta)*cos(theta) ...
            + (Qtheta+GGtheta)/(mstar*L^2*cos(phi)^2);
        ddphi = -((dtheta+dta)^2 ...
            + 3*(dta^2/(1+e*cos(ta)))*cos(theta)^2)*sin(phi)*cos(phi) ...
            + (Qphi+GGphi)/(mstar*L^2);

        dI = 0;
        
        dstate = [da; de; di; dRAAN; daop; dta; dtheta; dphi; ddtheta; ddphi; dI];
    end

    function [Bx, By, Bz] = MagField_NonTilted(states)
        % Earth’s magnetic field vector in the Euler–Hill frame whose 
        % components are defined by (assuming a nontilted dipole)
        a = states(1);
        e = states(2);
        i = states(3);
        RAAN = states(4);
        aop = states(5);
        ta = states(6);
        p = a*(1-e^2);
        
        Bo = 3.12e-5;
        Bo_over_R3 = Bo/(6378/(p/(1+e*cos(ta))))^3;

        Bx = -2*Bo_over_R3*sin(aop + ta)*sin(i);
        By = Bo_over_R3*cos(aop+ta)*sin(i);
        Bz = Bo_over_R3*cos(i);
    end

    function [fr, ftheta, fh] = edt_forces(states, tether_param)
        theta = states(7);
        phi = states(8);
        I = states(11);

        L = tether_param(1);
        m = sum(tether_param(2:4));

        [Bx, By, Bz] = MagField_NonTilted(states);

        fr = I*L/m*(Bz*sin(theta)*cos(phi)-By*sin(phi));
        ftheta = I*L/m*(Bx*sin(phi)-Bz*cos(theta)*cos(phi));
        fh = I*L/m*(By*cos(theta)*cos(phi)-Bx*sin(theta)*cos(phi));
    end

    function [Qtheta, Qphi] = edt_torque(states, tether_param)
        % assumes a uniform current flowing in the tether. In the case 
        % where the mass distribution is perfectly symmetrical, then the 
        % net torque about the center of mass is zero. This is an ideal 
        % scenario that is unlikely to be achieved in practice. This effect
        % will still be present in bare-wire tethers, but the relationship
        % is not as straightforward because of the nonuniform distribution 
        % of the electric current.
        theta = states(7);
        phi = states(8);
        I = states(11);

        L = tether_param(1);
        m1 = tether_param(2);
        m2 = tether_param(3);
        mt = tether_param(4);
        m = sum(tether_param(2:4));

        % nondimensional parameter that defines a measure of the average 
        % moment arm of the electromagnetic torque acting on the tether 
        % relative to the system center of mass
        PHI = (m1^2-m2^2+mt*(m1-m2))/m^2; 
        [Bx, By, Bz] = MagField_NonTilted(states);

        Qtheta = (I*L^2/2)*PHI*cos(phi)*...
            (sin(phi)*(Bx*cos(theta)+By*sin(theta))-Bz*cos(phi));
        Qphi = -(I*L^2/2)*PHI*(By*cos(theta)-Bx*sin(theta));
        Qtheta = 0;
        Qphi = 0;
    end

    function [GGtheta, GGphi] = gravity_grad_torque(states, tether_param, mu)
        a = states(1);
        e = states(2);
        i = states(3);
        RAAN = states(4);
        aop = states(5);
        ta = states(6);
        theta = states(7);
        phi = states(8);

        L = tether_param(1);
        m1 = tether_param(2);
        m2 = tether_param(3);
        mt = tether_param(4);
        m = sum(tether_param(2:4));

        [ r , v ] = coes2state( [sqrt(mu*a*(1-e^2)), i, e, RAAN, aop, ta] , mu );

        % need L in terms ECI to get r new
        GGtheta = 0;
        GGphi = 0;
%         grav_m1 = -mu*m1/r^3;
%         grav_m2 = -mu*m2/(r-L*cos)^3;
    end

    function [ r , v ] = coes2state( COES , mu )
    % h , inc , ecc , RAAN , omega , theta , a , rp , ra 
            
        h = COES(1) ;
        inc = COES(2) ;
        ecc = COES(3) ;
        RAAN = COES(4) ;
        omega = COES(5) ;
        theta = COES(6) ;
    
        r_peri = (h^2/mu) * ( 1/( 1 + ecc*cos(theta) ) ) * [ cos( theta ) ; sin( theta ) ; 0 ] ;
        v_peri = (mu/h) * [ -sin( theta ) ; ecc+cos(theta) ; 0 ] ;
    
        Q(1,1) = -sin(RAAN)*cos(inc)*sin(omega) + cos(RAAN)*cos(omega) ;
        Q(1,2) = -sin(RAAN)*cos(inc)*cos(omega) - cos(RAAN)*sin(omega) ;
        Q(1,3) = sin(RAAN)*sin(inc) ;
        Q(2,1) = cos(RAAN)*cos(inc)*sin(omega) + sin(RAAN)*cos(omega) ;
        Q(2,2) = cos(RAAN)*cos(inc)*cos(omega) - sin(RAAN)*sin(omega) ;
        Q(2,3) = -cos(RAAN)*sin(inc) ;
        Q(3,1) = sin(inc)*sin(omega) ;
        Q(3,2) = sin(inc)*cos(omega) ;
        Q(3,3) = cos(inc) ;
    
        r = Q*r_peri ;
        v = Q*v_peri ;
    end
end