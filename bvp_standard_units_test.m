clear; close all; clc;

m1 = 10; % lower mass
m2 = 10; % upper mass
mt = 10; % tether mass
L = 5000; % tether length

h0 = 1000;
solinit = bvpinit(linspace(0,L,100),@mat4init,h0);
sol = bvp4c(@mat4ode, @mat4bc, solinit);

fprintf('h0 is approximately %7.3f.\n',...
            sol.parameters)

hint = linspace(0,L);
Shint = deval(sol,hint);

plot(hint,Shint)
% axis([0 L -4 4])
xlabel('h')
ylabel('y')
legend('y','y''')


function dstate = mat4ode(h,state,h0) % equation being solved
    e = 1.60217663e-19; % electron charge
    ninf = 2.0208e+11; % ionospheric plasma density
    p = .0005*2*pi; %perimeter in meters
    me = 9.1093837e-31; % mass of electron
    sigma = 3.5e7; % from paper
    At = .0005^2*pi;
    Em = 165; % Induced electric field; dot(cross(v, B), u) 
    mu = 1/172; % apparently for atomic oxygen
    gamma = .15e-3;

    Ie = state(1);
    PHI = state(2);
    
    dstate = zeros(2,1);
    if h < h0
        dstate(1) = e*ninf*p/pi*sqrt(((2*e)/me)*PHI);
    else
        dstate(1) = -e*ninf*p/pi*sqrt(((2*e)/me)*abs(PHI))*mum*(1+gamma*abs(PHI));
    end
    dstate(2) = Ie/(sigma*At)-Em;
end
%-------------------------------------------
function res = mat4bc(ya,yb,h0) % boundary conditions
res = [ya(1)
       yb(1)
       ya(1)-1];
end
%-------------------------------------------
function yinit = mat4init(x) % initial guess function
yinit = [cos(4*x)
        -4*sin(4*x)];
end

function I = oml_sbet(h,p,L,At,PHI)
    me = 9.1093837e-31; % mass of electron
    mi = 2.6561e-23;
    Em = 165; % Induced electric field; dot(cross(v, B), u)
    e = 1.60217663e-19; % electron charge
    sigma = 3.5e7; % from paper
    ht = 2*At/p; % characteristic transversal length
    ninf = 2.0208e+11; % ionospheric plasma density 
    mu = 1/172; % apparently for atomic oxygen
    gamma = .15e-3;
    
    delta = .5; %gamma*Em*L
    Lstar = (me*Em)^(1/3)/(e*2^(7/3))*(3*pi*sigma*ht/ninf)^2/3;
    lt = L/Lstar;
    % asymptotic expressions
    i_B = mu*lt^(3/2)/2*(1+3/5*delta)*(1-mu^(2/3)*(3/2)*(1+delta)*(1+3/5*delta)^(-1/3)...
        -mu*(1/70)*lt^(3/2)*(81*delta^2+165*delta+70)/(5+3*delta));
    curlye_B = (2*i_B)^2/3*(1+4/15*i_B+53/360*i_B^2+3479/35640*i_B^3+3037/42768*i_B^4);
    
    curlye = h/Lstar;
    Isc = sigma*Em*At;
    I = Isc*ie(curlye);
    phie = PHI/(Em*Lstar);
    
    % NEED TO SOLVE BVP

end