function [m, phi, LAMBDAt, hG, Is] = params_2_sbet_params(m1, m2, mt, L)
    % m2 is upper mass
    % m1 is lower mass
    % mt is tether mass
    % L is tether length

% (m1, m2, mt) => (m, phi, LAMBDAt)
    m = m1 + m2 + mt; % System mass
    LAMBDAt = mt/m; % nondimensional tether mass
    
    cos2_phi = 1/m*(m1+.5*mt); 
    sin2_phi = 1/m*(m2+.5*mt);

    phi = atan2(sqrt(sin2_phi),sqrt(cos2_phi)); % measure of mass distribution

    hG = L*cos(phi)^2; % position of center of mass below m2
    Is = 1/12*m*L^2*(3*sin(2*phi)^2-2*LAMBDAt); % moment of inertia
end