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