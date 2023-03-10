function [ t , states] = BasicTether( tspan , sc_state0, tether_state0, tether_param, mu , tol )
% Uses Classical Orbital Elements a (semimajor axis), e (eccentricity), 
% i (inclination), RAAN (right ascension of ascending node), aop (argument
% of periapsis), ta (true anomaly)
%
% Considers electrodynamic force and gravity gradient
% will at some point maybe: atmospheric drag, spherical harmonics,
% third body of sun and moon, srp
% tether is treated as dumbbell model, rigid with lumped masses at the end


    opts = odeset('RelTol', tol, 'AbsTol', tol ) ;
    [ t , states ] = ode45(@gauss_variations, tspan , ...
        [ sc_state0 ; tether_state0], opts, tether_param, mu) ;


    function dstate = gauss_variations(t, states, tether_param, mu)

        % set orbit states with friendly names
        a = states(1);
        e = states(2);
        i = states(3);
        RAAN = states(4);
        aop = states(5);
        ta = states(6);

        % set attitude states with friendly names
        roll = states(7);
        pitch = states(8);
        yaw = states(9);
        droll = states(10);
        dpitch = states(11);
        dyaw = states(12);
        
        % set tether parameters with friendly names
        L = tether_param(1)*1000;
        m1 = tether_param(2);
        m2 = tether_param(3);
        mt = tether_param(4);
        m = sum(tether_param(2:4));
        Ix = tether_param(5);
        Iy = tether_param(6);
        Iz = tether_param(7);
        current_type = tether_param(8);
        current_val = tether_param(9);

        if current_type == 0
            I = current_val;
        elseif current_type == 1
            I = OML(tether_param);
        elseif current_type == 2
            % while energy is less than required to reach limit from 
            limit_libration = 35;
            [I, ~] = controlled_current(limit_libration, current_val, roll, pitch, droll, dpitch, m, L);
        end

        % Find the instantaneous Lorentz force, Lorentz torque, and gravity
        % gradient
        [Bx, By, Bz] = MagField_NonTilted(states);
%         [Bx, By, Bz] = MagField_igrf(states, t, mu);
%          r      theta     h
        [f_yaw, f_roll, f_pitch] = edt_forces(states, tether_param, I, Bx, By, Bz);
        Tq = edt_torque(states, tether_param, I, Bx, By, Bz);
        Tgg = gravity_grad_torque(states, tether_param, mu);

        % Calculate supplementary orbital elements
        p = a*(1-e^2);
        h = sqrt(mu*p);
        r = norm(states(1:3));
        n = sqrt(mu/a^3);
        u = aop + ta;

        % orbital elements change using Vallado p.636
        da = (2/(n*sqrt(1-e^2)))*(e*sin(ta)*f_yaw + (p/r)*f_roll);
        de = (sqrt(1-e^2)/(n*a))*(sin(ta)*f_yaw + ...
            (cos(ta)+(e+cos(ta))/(1+e*cos(ta)))*f_roll);
        di = (r*cos(u)/(n*a^2*sqrt(1-e^2)))*f_pitch;
        dRAAN = (r*sin(u)/(n*a^2*sqrt(1-e^2)*sin(i)))*f_pitch;
        daop = (sqrt(1-e^2)/(n*a*e))*(-cos(ta)*f_yaw + ...
            sin(ta)*(1+r/p)*f_roll) - (r*cot(i)*sin(u)/h)*f_yaw;
        dta = h/r^2 + (1/(e*h))*(p*cos(ta)*f_yaw - (p+r)*sin(ta)*f_roll);
        
        % attitude change from 16.1 of Spacecraft dynamics
        body_torque = Tgg+Tq;
        ddroll = (body_torque(1) - (Iz - Iy)*dpitch*dyaw)/Ix;
        ddpitch = (body_torque(2) - (Ix - Iz)*droll*dyaw)/Iy;        
        ddyaw = 0;
        dstate = [da; de; di; dRAAN; daop; dta; droll; dpitch; dyaw; ddroll; ddpitch; ddyaw];
    end

end