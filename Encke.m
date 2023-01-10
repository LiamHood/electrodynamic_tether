function [ t , r , v ] = Encke( dt , tspan , r0 , v0 , mu , forces , A , m )
    dr = zeros(3,1) ;
    eps = 0;
    f = 0;
    rp = r0 ;
    vp = v0;
    r = r0 ;
    v = v0 ;
    t = tspan(1) ;
    ii = 2 ;
    dv = zeros(3,1) ;
    while t(ii-1) < tspan(2) && norm(r(:,ii-1)) >= (6378 + 100 ) 
%        
        [ r(:,ii) , v(:,ii) ] = NewStateUV( r(:,ii-1) , v(:,ii-1) , dt , mu  ) ; 
        eps = dot( r(:,ii) , dr(:,ii-1) )/norm( r(:,ii) )^2 ;
        if eps ~= 0
            f = ( 1/eps )*( 1 - ( 1/ ( 1 - 2*eps )^1.5 ) ) ;
        else
            f = 0 ;
        end
        if contains( forces , "drag" )
            apd = DragAcceleration( r(:,ii) , v(:,ii) , A , m ) ;
        else
            apd = zeros(3,1) ;
        end
        if contains( forces , "gravity" )
            if contains( forces , "J2" )
                aph = J2accel( r ) ; 
            end
            if contains( forces , "J3" )
                aph = aph + J3accel( r ) ; 
            end
        else 
            aph = zeros(3,1) ; 
        end
        if contains( forces , "nbody" )
            apn = zeros(3,1) ; 
        else
            apn = zeros(3,1) ;
        end
        if contains( forces , "srp" )
            aps = zeros(3,1) ; 
        else
            aps = zeros(3,1) ;
        end
        ap(:,ii) = apd + aph + apn + aps ;
        da(:,ii) = ap(:,ii) + ( mu/norm(r(:,ii))^3 )*( f*eps*rp(:,ii-1) - dr(:,ii-1) ) ;
        dv = da(:,ii)*dt + dv ;
        t(end)/(60*60)
        dr(:,ii) = .5*da(:,ii)*dt^2 +  dv*dt + dr(:,ii-1) ;
        rp(:,ii) = r(:,ii) + dr(:,ii) ;
        vp(:,ii) = v(:,ii) + dv ;
        abs(norm(dr(ii))/norm(rp))*100
        if abs(norm(dr(ii))/norm(rp)) > 1e-4
            r(:,ii) = rp(:,ii) ;
            v(:,ii) = vp(:,ii) ;
            dr = zeros(3,1) ;
            dv = zeros(3,1) ;
        end  
        t(ii) = t(ii-1) + dt ;
        ii = ii + 1 ;
    end
    
    for ii = 1:length(t)
        mag_dr(ii) = norm(dr(:,ii));
        mag_ap(ii) = norm(ap(:,ii));
        t_mod(ii) = t(ii)/(60*60);
    end
    figure
    plot(t_mod, mag_dr)

    figure
    plot(t_mod, mag_dr)

    figure
    plot3(r(1,:),r(2,:), r(3,:))
    hold on
    plot3(rp(1,:), rp(2,:), rp(3,:))
    hold off

        function aP = J2accel( r )
            mu = 398600 ;
            J2 = -1.08262617385222e-3 ;
            Re = 6378 ;
            aP = zeros(3,1) ;
            rmag = norm(r) ;
            aP(1) = -3*J2*mu*Re^2*r(1)/(2*rmag^5)*(1-5*r(3)^2/rmag^2) ;
            aP(2) = -3*J2*mu*Re^2*r(2)/(2*rmag^5)*(1-5*r(3)^2/rmag^2) ;
            aP(3) = -3*J2*mu*Re^2*r(3)/(2*rmag^5)*(3-5*r(3)^2/rmag^2) ;
        end
    
        function aP = J3accel( r )
            mu = 398600 ;
            J3 = 2.53241051856772e-6 ;
            Re = 6378 ;
            aP = zeros(3,1) ;
            rmag = norm(r) ;
            aP(1) = -5*J3*mu*Re^3*r(1)/(2*rmag^7)*(3*r(3)-7*r(3)^3/rmag^2) ;
            aP(2) = -5*J3*mu*Re^3*r(2)/(2*rmag^7)*(3*r(3)-7*r(3)^3/rmag^2) ;
            aP(3) = -5*J3*mu*Re^3/(2*rmag^7)*(6*r(3)^2-7*r(3)^4/rmag^2 - (3/5)*rmag^2) ;
        end
    
        

end