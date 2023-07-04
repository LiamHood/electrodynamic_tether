clear; close all; clc;

n = 1e2;
mm = 1e2;

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
m_additional = 100;
% figure
% plot(L*1e-3, rad2deg(phi))
% xlabel('Length [km]')
% ylabel('PHI [degree]')

ninf = linspace(4e2,1.8e6,mm).*1e6;
% ninf = logspace(2, 6,mm).*1e6;

for ii = 1:n
    for jj = 1:mm      
        mt(ii) = mt_per_m*L(ii);
        mt_coiled = mt_tot - mt(ii) ;
        m1(ii) = m_low ;
        m2(ii) = m_up + m_additional + mt_coiled;
        [m, phi(ii), LAMBDAt(ii), ~, ~] = params_2_sbet_params(m1(ii), m2(ii), mt(ii), L(ii));
        % Varies on one orbit
%         ninf = 2.0208e+11; % ionospheric plasma density 
        Em = 165e-3; % V/m
        delta = gamma1*Em*L(ii);
        Lstar = (((me*Em)^(1/3))/(e*2^(7/3)))*(3*pi*sigma*ht/ninf(jj))^(2/3);
        lt = L(ii)/Lstar;
    
        % non-dimensional parameters
        eps0 = Em/L(ii)*(12*LAMBDAt(ii)/(3*sin(2*phi(ii))^2-2*LAMBDAt(ii)))*(mum/mue_si)*(sigma/rhov);
        fhat_approx = mu*lt^1.5*.3*(1+5/7*delta)*(cos(phi(ii))^2-(5/18)*((9+7*delta)/(7+5*delta)));
        epsilon(ii,jj) = eps0*fhat_approx;

        phi_star(ii) = acosd(sqrt((5/18)*((9+7*delta)/(7+5*delta))));
    end
end

figure
hold on
plot(L*1e-3, rad2deg(phi))
plot(L*1e-3, phi_star)
hold off
xlabel('Length of Tether [km]')
ylabel('\Phi [degree]')
legend('\Phi actual', '\Phi^*')

figure
plot(L*1e-3, LAMBDAt)
xlabel('Length of Tether [km]')
ylabel('\Lambda_t')

figure
surf(ninf, L*1e-3, epsilon, abs(epsilon), 'EdgeColor', 'none')
xlabel('N_infinite [1/m^3]')
ylabel('Length of Tether [km]')
zlabel('\epsilon')


figure
surf(ninf, rad2deg(phi), epsilon, abs(epsilon), 'EdgeColor', 'none')
xlabel('N_infinite [1/m^3]')
ylabel('\Phi [degree]')
zlabel('\epsilon')

figure
hold on
plot(L*1e-3, epsilon(:,1))
plot(L*1e-3, epsilon(:,ceil(mm/4)))
plot(L*1e-3, epsilon(:,2*ceil(mm/4)))
plot(L*1e-3, epsilon(:,end))
legend('electron density = 4e2 m^-^2', 'electron density = 4e5 m^-^2', ...
    'electron density = 1e6 m^-^2', 'electron density = 1.8e6 m^-^2')
xlabel('Length of Tether [km]')
ylabel('\epsilon')

figure
hold on
plot(rad2deg(phi), epsilon(:,1))
plot(rad2deg(phi), epsilon(:,ceil(mm/4)))
plot(rad2deg(phi), epsilon(:,2*ceil(mm/4)))
plot(rad2deg(phi), epsilon(:,end))
legend('electron density = 4e2 m^-^2', 'electron density = 4e5 m^-^2', ...
    'electron density = 1e6 m^-^2', 'electron density = 1.8e6 m^-^2')
xlabel('\Phi [degree]')
ylabel('\epsilon')