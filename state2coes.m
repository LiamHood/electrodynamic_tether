function COES = state2coes( state, mu )
% all angles output in radians
% h , inc , ecc , RAAN , omega , theta , a , rp , ra 
    
    Kh = [ 0 0 1 ] ; % K hat
    R = state(1:3) ;
    V = state(4:6) ;
    r = norm( R ) ;
    v = norm( V ) ;
    vr = dot( R , V )/r ; % radial velocity
    H = cross( R , V ) ; % specific angular momentum
    h = norm( H ) ; % specific angular momentum
    inc = acos( H(3)/h ) ; %inclination
    ECC = (1/mu)*( ( v^2 - mu/r )*R - r*vr*V ) ; %eccentricity vector
    ecc = norm( ECC ) ; % eccentricity
    N = cross( Kh , H ) ; % Node line
    n = norm( N ) ;
if n ~= 0
    RAAN = acos(N(1)/n) ; %Right ascension of ascending node
    if N(2) < 0 
        RAAN = 2*pi - RAAN ; %Right ascension of ascending node
    end
else
    RAAN = 0 ;
end 
    
if n ~= 0 
    if ecc >= 0 
        omega = acos(dot(N,ECC)/(n*ecc)) ; % Argument of perigee
        if ECC(3) < 0 
            omega = 2*pi - omega ; % Argument of perigee
        end
    else
        omega = 0 ;
    end
else
    omega = 0 ;
end
    
if ecc > 0
    theta = acos( dot( ECC , R )/( ecc*r ) ) ;         
    if vr < 0 
        theta = 2*pi - theta ; 
    end
else
    cp = cross( N , R ) ;
    if cp(3) >= 0
        theta = acos( dot( N , R )/( n*r ) ) ;
    else
        theta = 2*pi - acos( dot( N , R )/( n*r ) ) ;
    end
end

    a = (h^2)/( mu*( 1 - ecc^2 ) ) ; % semi-major axis
    rp = a*( 1 - ecc ) ;
    ra = a*( 1 + ecc ) ;

    
    COES = [ h , inc , ecc , RAAN , omega , theta , a , rp , ra ] ;
end


