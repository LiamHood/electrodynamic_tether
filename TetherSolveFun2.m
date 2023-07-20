function F = TetherSolveFun2(s0, tether, constants, parameters) % boundary conditions  
%     p = tether.perimeter; %perimeter in meters
%     At = tether.At;
%     me = constants.me; % mass of electron
%     Em = ; % Induced electric field; dot(cross(v, B), u)
%     e = parameters(5); % electron charge
%     sigma = parameters(6); % from paper
%     ninf = parameters(7); % ionospheric plasma density 
%     mu = parameters(8); % apparently for atomic oxygen
%     gamma = parameters(9);

    Ie_0 = 0;
    PHIe_0 = s0(1); 
    h0 = s0(2);
    state0 = [Ie_0, PHIe_0];
    h_span1 = [0, h0] ;
    h_span2 = [h0, tether.L] ;
    opts = odeset( 'RelTol' , 1e-8 , 'AbsTol' , 1e-8 ) ;
    [t1, s1] = ode45( @bvp_ode_low, h_span1, state0, opts, tether, constants, parameters) ;
    sB = s1(end,:);
    [t2, s2] = ode45( @bvp_ode_high, h_span2, sB, opts, tether, constants, parameters) ;
    sf = s2(end,:);
    % [ie, phie]

    F(1) = sB(1)-max(max(s1(:,1)), max(s2(:,1))); % Ie = max current = I_B, at B
    F(2) = sB(2); % PHIe = 0 at B
    F(3) = sf(1)*100; % Ie = 0, I_B, at B
%     F(4) = sum(abs(s1(:,2))-s1(:,2)); % bias should be positive until h_0
    
%     fprintf("Ie max at middle: %f\n",F(1))
%     fprintf("PHIe zero at middle: %f\n",F(2))
%     fprintf("Ie zero at end: %f\n",F(3))



    function dstate = bvp_ode_low(h, state, tether, constants, parameters) % equation being solved
        Ie = state(1);
        PHIe = state(2);
          
        dstate = zeros(2,1);
        dstate(1) = constants.e*parameters.ninf*tether.perimeter/pi*sqrt((2*constants.e/constants.me)*abs(PHIe));
        dstate(2) = Ie/(tether.sigma*tether.At)-parameters.Em;
    end

    function dstate = bvp_ode_high(h, state, tether, constants, parameters) % equation being solved
        Ie = state(1);
        PHIe = state(2);
    
        dstate = zeros(2,1);
        dstate(1) = -constants.e*parameters.ninf*tether.perimeter/pi*sqrt((2*constants.e/constants.me)*abs(PHIe))*constants.mu*(1+constants.gamma*abs(PHIe));
        dstate(2) = Ie/(tether.sigma*tether.At)-parameters.Em;
    end
end