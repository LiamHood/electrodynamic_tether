function B = sbet_magfield(a, i, ta)
    B = zeros(3,1);
    mu = 8e15; % Tesla*m^3
    a = a*1e3; % change to meters from km
    B(1) = -2*mu/a^3*sin(i)*sin(ta);
    B(2) = -mu/a^3*cos(i);
    B(3) = mu/a^3*sin(i)*cos(ta);
end