function F = TetherSolveFun(s0, parameters, L) % boundary conditions  
    p = parameters(1); %perimeter in meters
    At = parameters(2);
    me = parameters(3); % mass of electron
    Em = parameters(4); % Induced electric field; dot(cross(v, B), u)
    e = parameters(5); % electron charge
    sigma = parameters(6); % from paper
    ninf = parameters(7); % ionospheric plasma density 
    mu = parameters(8); % apparently for atomic oxygen
    gamma = parameters(9);

    Ie_0 = 0;
    PHIe_0 = s0(1); 
    I_B = s0(2);
    h0 = s0(3);
    state0 = [Ie_0, PHIe_0];
    h_span1 = [0, h0] ;
    h_span2 = [h0, L] ;
    opts = odeset( 'RelTol' , 1e-8 , 'AbsTol' , 1e-8 ) ;
    [t1, s1] = ode45( @bvp_ode_low, h_span1, state0, opts) ;
    sB = s1(end,:);
    [t2, s2] = ode45( @bvp_ode_high, h_span2, sB, opts) ;
    sf = s2(end,:);
    % [ie, phie]

    F(1) = sB(1)-I_B; % Ie = max current = I_B, at B
    F(2) = sB(2); % PHIe = 0 at B
    F(3) = sf(1)*10; % Ie = 0, I_B, at B
    F(4) = imag(s0(1));
    F(5) = imag(s0(2));
    F(6) = imag(s0(3));

    function dstate = bvp_ode_low(h,state) % equation being solved
        Ie = state(1);
        PHIe = state(2);
          
        dstate = zeros(2,1);
        dstate(1) = e*ninf*p/pi*sqrt((2*e/me)*PHIe);
        dstate(2) = Ie/(sigma*At)-Em;
    end

    function dstate = bvp_ode_high(h,state) % equation being solved
        Ie = state(1);
        PHIe = state(2);
    
        dstate = zeros(2,1);
        dstate(1) = -e*ninf*p/pi*sqrt((2*e/me)*abs(PHIe))*mu*(1+gamma*abs(PHIe));
        dstate(2) = Ie/(sigma*At)-Em;
    end
end