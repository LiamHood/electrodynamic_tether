function I = OML(tether_param)
    % Finds the OML current
    e = 1.60217663e-19 ; % charge of an electron
    me = 9.1093837e-31; % mass of an electron
    R = tether_param(10); % Tether radius (m)
    L = tether_param(1)*1000; % Tether Length (m)
    p = R*2*pi; % Tether cross sectional parameter (m)
    phiP = 100; % Cylindrical probe bias
    N0 = 0; % ambient electron density
    
    
    I = e*N0*(L*p/pi)*sqrt(2*e*phiP/me);


end
    