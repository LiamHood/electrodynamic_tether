function [ r , v ] = NewStateUV( r0 , v0 , dt , mu )
% Find new position and velocity at some time in orbit by lagrange
% variables

    % Initial values
    r0mag = norm( r0 ) ; % radius in km
    v0mag = norm( v0 ) ; % speed in km/s
%     vr0 = dot( r0 , v0 ) / r0mag ; % radial speed (km/s)   
    visviva = v0mag^2/2 - mu/r0mag ;
    alpha = (2/r0mag) - ( v0mag^2 / mu ) ; % reciporical of semimajor axis
    if alpha > 1e-6 
        chi0 = sqrt( mu ) * dt * alpha ; % Initial universal anomaly guess
    elseif alpha < -1e-6
        a = 1/alpha ;
        chi0 = sign( dt )*sqrt( -a )*log( ( -2*mu*alpha*dt )/( dot( r0 , v0 ) + sign( dt )*sqrt( -mu*a )*( 1 - r0mag*alpha ) ) ) ;
    else 
        h = cross( r0 , v0 ) ;
        p = norm( h )^2 / mu ;
        s = 2*acot( 3*sqrt( mu / p^3 )*dt ) ;
        w = nthroot( atan( tan( s ) ) , 3 ) ;
        chi0 = sqrt( p )*2*cot( 2*w ) ;
    end
    
        % Newtons for UV
        ii = 1 ;
        err = 1 ;
        tol = 1e-8 ;
        lim = 1e4 ;
            while err >= tol
                z = alpha*chi0^2 ;    
                [ c3 , c2 ] = StumpfCalc( z ) ;
                bottom = chi0^2*c2 + ( dot( r0 , v0 )/sqrt( mu ) )*chi0*( 1 - z*c3 ) + r0mag*( 1 - z*c2 ) ;
                chi = chi0 + ( sqrt( mu )*dt - chi0^3*c3 - ( dot( r0 , v0 )/sqrt( mu ) )*chi0^2*c2 - r0mag*chi0*( 1 - z*c3 ) )/bottom ;
                err = abs( chi - chi0 ) ;
                chi0 = chi ;
                    if ii > lim
                        error([ 'Ran ' , num2str( lim ) , ' times without a solution' ])
                    end
            end   
    
    % Lagrange variables
    
    % new r
    f = 1 - ( chi^2/r0mag )*c2 ;
    g = dt - ( chi^3/sqrt( mu ) )*c3 ;
    r = f*r0 + g*v0 ; % new position
    rmag = norm( r ) ; % radius
    
    % new v
    fdot = ( sqrt(mu) / ( r0mag*rmag ) )*chi*( z*c3 - 1 ) ; %(alpha*chi^3*S-chi)
    gdot = 1 - ( chi^2/rmag )*c2 ; 
    v = fdot*r0 + gdot*v0 ; % new velocity
      
end