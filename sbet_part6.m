clear; close all; clc;

n = 1e2;
mm = 1e3;

% Constants
sigma = 3.5e7; % 1/(Ohm*m)
rhov = 2700; % kg/m^3
mue = 398600; 
mue_si = mue*1e9; % mu of earth in meters
mum = 8e15; % magnetic constant of earths mag field
me = 9.1093837e-31; % mass of electron
mi = 2.6561e-26; % mass of atomic oxygen
e = 1.60217663e-19; % electron charge
mu = sqrt(me/mi); % Ratio (me/mi)^1/2
At = 2.1096e-6;
p = 24e-3;
ht = 2*At/p;
gamma1 = .15e-3;

Lmax = 20e3;
L = linspace(1e2, Lmax, n); % m
m_max= 2e2 ; % kg
mt_tot = 114.9;
mt_per_m = (mt_tot/Lmax);
m_low = 100;
m_up = 100;
m_additional = linspace(0, m_max, mm);
for ii = 1:n
    for jj = 1:mm
        mt(ii,jj) = mt_per_m*L(ii);
        mt_coiled = mt_tot - mt(ii,jj) ;
        m1(ii,jj) = m_low ;
        m2(ii,jj) = m_up + m_additional(jj) + mt_coiled;
        [m, phi(ii,jj), LAMBDAt(ii,jj), ~, ~] = params_2_sbet_params(m1(ii,jj), m2(ii,jj), mt(ii,jj), L(ii));
    end
end
% figure
% plot(L*1e-3, rad2deg(phi))
% xlabel('Length [km]')
% ylabel('PHI [degree]')
figure
surf(m_additional, L*1e-3, rad2deg(phi), abs(rad2deg(phi)), 'EdgeColor', 'none')
xlabel('Mass of Sat to Deorbit [kg]')
ylabel('Length of Tether [km]')
zlabel('\Phi [degree]')

figure
surf(m_additional, L*1e-3, LAMBDAt, abs(LAMBDAt), 'EdgeColor', 'none')
xlabel('Mass of Sat to Deorbit [kg]')
ylabel('Length of Tether [km]')
zlabel('\Lambda_t')


for ii = 1:n
    for jj = 1:mm
%     phid = 10; % degrees
%     phi = deg2rad(phid);
%     LAMBDAt = mt/m;
%     [m1, m2, mt] = sbet_params_2_masses(m, phi, LAMBDAt);
            
        % Need LAMBDAt, phi
    
        
        
        % Varies on one orbit
        ninf = 2.0208e+11; % ionospheric plasma density 
        Em = 165e-3; % V/m
        delta = gamma1*Em*L(ii);
        Lstar = (((me*Em)^(1/3))/(e*2^(7/3)))*(3*pi*sigma*ht/ninf)^(2/3);
        lt = L(ii)/Lstar;
    
        % non-dimensional parameters
        eps0 = Em/L(ii)*(12*LAMBDAt(ii,jj)/(3*sin(2*phi(ii,jj))^2-2*LAMBDAt(ii,jj)))*(mum/mue_si)*(sigma/rhov);
        fhat_approx = mu*lt^1.5*.3*(1+5/7*delta)*(cos(phi(ii,jj))^2-(5/18)*((9+7*delta)/(7+5*delta)));
        epsilon(ii,jj) = eps0*fhat_approx;
    end
end

figure
surf(m_additional, L*1e-3, epsilon, abs(epsilon), 'EdgeColor', 'none')
xlabel('Mass of Sat to Deorbit [kg]')
ylabel('Length of Tether [km]')
zlabel('\epsilon')


figure
hold on
plot(L*1e-3, epsilon(:,1))
plot(L*1e-3, epsilon(:,ceil(mm/4)))
plot(L*1e-3, epsilon(:,2*ceil(mm/4)))
plot(L*1e-3, epsilon(:,3*ceil(mm/4)))
plot(L*1e-3, epsilon(:,end))
plot([L(1),L(end)]*1e-3, [.0005, .0005],'r')
plot([L(1),L(end)]*1e-3, [-.0005, -.0005], 'r')
legend('deorbit mass = 0', 'deorbit mass = 50', 'deorbit mass = 100', ...
    'deorbit mass = 150', 'deorbit mass = 200', 'Location', 'northwest')

