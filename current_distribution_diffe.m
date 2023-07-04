function dstate = current_distribution_diffe(h,state,parameters) 
    Ie = real(state(1));
    PHIe = real(state(2));

    p = parameters(1); %perimeter in meters
    At = parameters(2);
    me = parameters(3); % mass of electron
    Em = parameters(4); % Induced electric field; dot(cross(v, B), u)
    e = parameters(5); % electron charge
    sigma = parameters(6); % from paper
    ninf = parameters(7); % ionospheric plasma density 
    mu = parameters(8); % apparently for atomic oxygen
    gamma = parameters(9);
    h0 = parameters(10);

    dstate = zeros(2,1);
    if real(h) < h0
        dstate(1) = e*ninf*p/pi*sqrt((2*e/me)*abs(PHIe));
    else
        dstate(1) = -e*ninf*p/pi*sqrt((2*e/me)*abs(PHIe))*mu*(1+gamma*abs(PHIe));
    end
    dstate(2) = Ie/(sigma*At)-Em;
end