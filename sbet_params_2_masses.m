function [m1, m2, mt] = sbet_params_2_masses(m, phi, LAMBDAt)
% m: total mass;
% phi: mass distribution;
% LAMBDAt: tether mass/total mass
    m1 = m*(cos(phi)^2-.5*LAMBDAt);
    m2 = m*(sin(phi)^2-.5*LAMBDAt);
    mt = m - m1 - m2;

end