function dstate = sbet_libration(t, state, mu, epsilon)
    a = state(1);
    e = state(2);
    i = state(3);
    RAAN = state(4);
    aop = state(5);
    ta = state(6);
    theta = state(7); % in plane
    phi = state(8); % out of plane
    dtheta = state(9);
    dphi = state(10);

    
     % Calculate supplementary orbital elements
    p = a*(1-e^2);
    h = sqrt(mu*p);
    r = p/(1+e*cos(ta));
    n = sqrt(mu/a^3);
    u = aop + ta;

    f_roll = 0;
    f_pitch = 0;
    f_yaw = 0;
    da = (2/(n*sqrt(1-e^2)))*(e*sin(ta)*f_yaw + (p/r)*f_roll);
    de = (sqrt(1-e^2)/(n*a))*(sin(ta)*f_yaw + ...
        (cos(ta)+(e+cos(ta))/(1+e*cos(ta)))*f_roll);
    di = (r*cos(u)/(n*a^2*sqrt(1-e^2)))*f_pitch;
    dRAAN = (r*sin(u)/(n*a^2*sqrt(1-e^2)*sin(i)))*f_pitch;
    daop = (sqrt(1-e^2)/(n*a*e))*(-cos(ta)*f_yaw + ...
        sin(ta)*(1+r/p)*f_roll) - (r*cot(i)*sin(u)/h)*f_yaw;
    dta = h/r^2 + (1/(e*h))*(p*cos(ta)*f_yaw - (p+r)*sin(ta)*f_roll);
    ddtheta = 2*(1+dtheta)*dphi*tan(phi) - 1.5*sin(2*theta)...
        -epsilon*(sin(i)*tan(phi)*(2*sin(ta)*cos(theta) - cos(ta)*sin(theta)) + cos(i));
    ddphi = -sin(phi)*cos(phi)*((1+dtheta)^2 + 3*cos(theta)^2)...
        +epsilon*sin(i)*(2*sin(ta)*sin(theta)+cos(ta)*cos(theta));
    
    dstate = [da; de; di; dRAAN; daop; dta; dtheta; dphi; ddtheta; ddphi];
end