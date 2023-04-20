function [ t , states] = BasicTether( tspan , sc_state0, tether_state0, tether_param, mag_model, mu , tol , net)
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
        [ sc_state0 ; tether_state0], opts, tether_param, mag_model, mu, net) ;


    function dstate = gauss_variations(t, states, tether_param, mag_model, mu, net)
        jdate = 2458849.5 + t/(24*3600);
        t/3600
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

        if mag_model == 1
            [Bx, By, Bz] = MagField_igrf(states, t, mu);
        else
            [Bx, By, Bz] = MagField_NonTilted(states); % a, e, i, RAAN, omega, theta
        end

        if current_type == 0
            I = current_val;
        elseif current_type == 1
            % while energy is less than required to reach limit from 
            limit_libration = 35;
            [I, ~] = controlled_current(limit_libration, current_val, roll, pitch, droll, dpitch, m, L);
        elseif current_type == 2
            N0 = 2.0208e5;
            I = OML(tether_param, N0*(1e2)^3);
        elseif current_type == 3
            utc_time = datetime(jdate,'ConvertFrom','juliandate');
            [ r , ~ ] = coes2state([sqrt(mu*a*(1-e^2)), i, e, RAAN, aop, ta], mu );
            lla = eci2lla(r'*1e3,[year(utc_time), month(utc_time), day(utc_time),...
                hour(utc_time), minute(utc_time),second(utc_time)]);
            [zd,~,~] = timezone(lla(2));
            time = zd + hour(utc_time) + minute(utc_time)/60 + second(utc_time)/(60*60);
            if time < 0
                time = time + 24;
            elseif time > 24
                time = time - 24;
            end

            input_raw = [jdate - 2433282.5; time; lla(1); lla(2); lla(3)*1e-3];
            input_avg = [25581.9999999703, 10.3747458975316, 0.284981889543736, 180.007038172997, 1049.08650750044];
            input_std = [8.66045976272869, 6.02542658129692, 52.0528796594471, 103.910546715735, 548.969414730606];
            for ii = 1:length(input_raw)
                input(ii) = (input_raw(ii)-input_avg(ii))/input_std(ii);
            end
            N0_raw = net(input');
            N0_std = 1.522683511834375e5;
            N0_avg = 8.668429073888397e4;
            N0 = abs(N0_raw*N0_std+N0_avg);
            I = OML(tether_param, N0*(1e2)^3);
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
        r = p/(1+e*cos(ta));
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