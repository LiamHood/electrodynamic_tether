function [ t , r , v ] = TwoBody( tspan , r0 , v0 , mu , tol )
    opts = odeset( 'RelTol' , tol , 'AbsTol' , tol ) ;
    [ t , rv ] = ode45( @TwoBodyForceFun , tspan , [ r0 ; v0 ] , opts , mu ) ;
    r = rv(:,1:3)' ;
    v = rv(:,4:6)' ;
    
    function drdv = TwoBodyForceFun( t , rv , mu )
        rvec = rv(1:3) ;
        vvec = rv(4:6) ;

        dr = vvec ;
        dv = ( -mu / norm(rvec)^3 )*rvec ;
        drdv = [ dr ; dv ] ;
    end

end