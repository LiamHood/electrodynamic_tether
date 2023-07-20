function dstate = current_distribution_diffe2(h,state, h0, tether, constants, parameters) 
    Ie = state(1);
    PHIe = state(2);

    dstate = zeros(2,1);
    if h < h0
        dstate(1) = constants.e*parameters.ninf*tether.perimeter/pi*sqrt((2*constants.e/constants.me)*abs(PHIe));
    else
        dstate(1) = -constants.e*parameters.ninf*tether.perimeter/pi*sqrt((2*constants.e/constants.me)*abs(PHIe))*constants.mu*(1+constants.gamma*abs(PHIe));
    end
    dstate(2) = Ie/(tether.sigma*tether.At)-parameters.Em;
end