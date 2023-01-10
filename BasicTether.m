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
        phi = states(7);
        theta = states(8);
        psi = states(9);
        dphi = states(10);
        dtheta = states(11);
        dpsi = states(12);
        
        % set tether parameters with friendly names
        L = tether_param(1);
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
            I = current_val;
            lim_e = m*L*1000*(1-cos(deg2rad(limit_libration)));
            pe_p = m*L*1000*(1-cos(phi));
            ke_p = (1/2)*Ix*(dphi*4)^2;
%             ke_p = (1/2)*Ix*(dphi)^2;
            if phi < 0
                pe_p = -pe_p;
            end
            if dphi < 0
                ke_p = -ke_p;
            end
            e_p = pe_p + ke_p;
            pe_t = m*L*1000*(1-cos(theta));
            ke_t = (1/2)*Iy*(dtheta*4)^2;
%             ke_t = (1/2)*Iy*(dtheta)^2;
            if theta < 0
                pe_t = -pe_t;
            end
            if dtheta < 0
                ke_t = -ke_t;
            end
            e_t = pe_t + ke_t;
            if abs(e_p) > lim_e
                I = 0;
            elseif abs(e_t) > lim_e
                I = 0;
            end
        end

        % Find the instantaneous Lorentz force, Lorentz torque, and gravity
        % gradient
%         [Bx, By, Bz] = MagField_NonTilted(states);
        [Bx, By, Bz] = MagField_igrf(states, t, mu);
        [fr, fs, fw] = edt_forces(states, tether_param, I, Bx, By, Bz);
        Tq = edt_torque(states, tether_param, I, Bx, By, Bz);
        Tgg = gravity_grad_torque(states, tether_param, mu);

        % Calculate supplementary orbital elements
        p = a*(1-e^2);
        h = sqrt(mu*p);
        r = norm(states(1:3));
        n = sqrt(mu/a^3);
        u = aop + ta;

        % orbital elements change using Vallado p.636
        da = (2/(n*sqrt(1-e^2)))*(e*sin(ta)*fr + (p/r)*fs);
        de = (sqrt(1-e^2)/(n*a))*(sin(ta)*fr + ...
            (cos(ta)+(e+cos(ta))/(1+e*cos(ta)))*fs);
        di = (r*cos(u)/(n*a^2*sqrt(1-e^2)))*fw;
        dRAAN = (r*sin(u)/(n*a^2*sqrt(1-e^2)*sin(i)))*fw;
        daop = (sqrt(1-e^2)/(n*a*e))*(-cos(ta)*fr + ...
            sin(ta)*(1+r/p)*fs) - (r*cot(i)*sin(u)/h)*fw;
        dta = h/r^2 + (1/(e*h))*(p*cos(ta)*fr - (p+r)*sin(ta)*fs);
        
        % attitude change from 16.1 of Spacecraft dynamics
        wo = -sqrt(mu/r^3);
        body_torque = Tgg+Tq;
        ddphi = dpsi*wo + ((Iz-Iy)*(wo^2*phi+wo*dpsi)+body_torque(1))/Ix;
        ddtheta = body_torque(2)/Iy;        
        ddpsi = 0;
        dstate = [da; de; di; dRAAN; daop; dta; dphi; dtheta; dpsi; ddphi; ddtheta; ddpsi];
    end

end