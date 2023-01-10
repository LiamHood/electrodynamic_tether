function [ t , states] = encke_tether( tspan , sc_state0, tether_state0, tether_param, mag_model, mu , tol )
% Uses enckes method of perturbation propagation,
%
% Considers electrodynamic force and gravity gradient
% will at some point maybe: atmospheric drag, spherical harmonics,
% third body of sun and moon, srp
% tether is treated as dumbbell model, rigid with lumped masses at the end

    dt = 10;
    
    % set starting values
    rp = sc_state0(1:3) ;
    r(:,1) = sc_state0(1:3) ;
    v(:,1) = sc_state0(4:6) ;
    t(1) = tspan(1) ;
    ii = 2 ;
    dr = zeros(3,1) ;
    dv = zeros(3,1) ;
    
    % set attitude states with friendly names
    phi(1) = tether_state0(1);
    theta(1) = tether_state0(2);
    psi(1) = tether_state0(3);
    dphi(1) = tether_state0(4);
    dtheta(1) = tether_state0(5);
    dpsi(1) = tether_state0(6);
    
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
    
    
    % loop until time ends or satellite gets to 100 km of earth
    while t(ii-1) < tspan(2) && norm(r(:,ii-1)) >= (6378 + 100 ) 
        % propagate unperturbed orbit
        [ r(:,ii) , v(:,ii) ] = NewStateUV( r(:,ii-1) , v(:,ii-1) , dt , mu  ) ; 
        
        % find eps and f for encke
        eps = dot( r(:,ii) , dr )/norm( r(:,ii) )^2 ;
        if eps ~= 0
            f = ( 1/eps )*( 1 - ( 1/ ( 1 - 2*eps )^1.5 ) ) ;
        else
            f = 0 ;
        end

        if current_type == 0 % constant
            I = current_val;
        elseif current_type == 1 % OML, still needs electron density
            I = OML(tether_param);
        elseif current_type == 2 % controlled
            limit_libration = 35; % limit degrees for calculating energy limit
            I = current_val;
            lim_e = m*L*1000*(1-cos(deg2rad(limit_libration)));

            % find "energy"
            pe_p = m*L*1000*(1-cos(phi(ii)));
            ke_p = (1/2)*Ix*(dphi(ii)*4)^2;
            if phi < 0
                pe_p = -pe_p;
            end
            if dphi < 0
                ke_p = -ke_p;
            end
            e_p = pe_p + ke_p;
            pe_t = m*L*1000*(1-cos(theta));
            ke_t = (1/2)*Iy*(dtheta*4)^2;
            if theta < 0
                pe_t = -pe_t;
            end
            if dtheta < 0
                ke_t = -ke_t;
            end
            e_t = pe_t + ke_t;

            % set to 0 if too much energy
            if abs(e_p) > lim_e
                I = 0;
            elseif abs(e_t) > lim_e
                I = 0;
            end
        end
    
        % Find COES so that the function developed for VOP work
        COES = state2coes([r(:,ii-1); v(:,ii-1)], mu ); % h , inc , ecc , RAAN , omega , theta , a , rp , ra 
        states_wCOES = [COES(7), COES(3), COES(2), COES(4), COES(5), COES(6), phi(ii-1), theta(ii-1), psi(ii-1), dphi(ii-1), dtheta(ii-1), dphi(ii-1)];
        
        % Calculate the magnetic field
        if mag_model == 1
            [Bx, By, Bz] = MagField_igrf(states_wCOES, t(ii-1), mu);
        else
            [Bx, By, Bz] = MagField_NonTilted(states_wCOES); % a, e, i, RAAN, omega, theta
        end
        % calculate lorentz acceleration, lorentz torque, and gravity gradient
        [fr, fs, fw] = edt_forces(states_wCOES, tether_param, I, Bx, By, Bz); % actually acceleration not force
        Tq = edt_torque(states_wCOES, tether_param, I, Bx, By, Bz);
        Tgg = gravity_grad_torque(states_wCOES, tether_param, mu);

        % convert to rsw
        r_rsw = r(:,ii-1)/norm(r(:,ii-1));
        w_rsw = cross(r(:,ii-1),v(:,ii-1))/norm(cross(r(:,ii-1),v(:,ii-1)));
        s_rsw = cross(w_rsw, r_rsw)/norm(cross(w_rsw, r_rsw));
        rsw2eci = [r_rsw, s_rsw, w_rsw];

        % calculate change in
        ap_rsw = [fr, fs, fw]'; % rsw acceleration
        ap = rsw2eci*ap_rsw; % eci acceleration
        da = ap + ( mu/norm(r(:,ii))^3 )*( f*eps*rp - dr ) ; % difference in acceleration
        dv = da*dt + dv ; % diff velocity
        dr = .5*da*dt^2 + dv*dt + dr ; % diff position
        rp = r(:,ii) + dr ; % position perturbed
        vp = v(:,ii) + dv ; % velocity perturbed
%         track_dr(ii) = norm(dr);        
        if norm(dr)/norm(rp) > 1e-4
            r(:,ii) = rp ;
            v(:,ii) = vp ;
            dr = zeros(3,1);
            dv = zeros(3,1);
        end  
        % attitude change from 16.1 of Spacecraft dynamics
        wo = -sqrt(mu/norm(r(:,ii))^3);
        body_torque = Tgg+Tq;
        ddphi = dpsi(ii-1)*wo + ((Iz-Iy)*(wo^2*phi(ii-1)+wo*dpsi(ii-1))+body_torque(1))/Ix;
        ddtheta = body_torque(2)/Iy;        
        ddpsi = 0;
        dphi(ii) = ddphi*dt + dphi(ii-1);
        dtheta(ii) = ddtheta*dt + dtheta(ii-1);
        dpsi(ii) = 0;
        phi(ii) = ddphi*dt^2 + dphi(ii-1)*dt + phi(ii-1);
        theta(ii) = ddtheta*dt^2 + dtheta(ii-1)*dt + theta(ii-1);
        psi(ii) = 0;

        t(ii) = t(ii-1) + dt ;
        ii = ii + 1 ;
    end
    states = [r; v; phi; theta; psi; dphi; dtheta; dpsi];
end