clear; close all; clc;
% https://arc-aiaa-org.proxy.wichita.edu/doi/pdf/10.2514/6.2004-5309

%% Define tether
tether.m1 = 321.6; % lower mass [kg]
tether.m2 = 563.5; % upper mass [kg]
tether.L = 20000; % tether length [km]
tether.density = 2700; % [kg/m^3]
tether.thickness = .18e-3; % [m]
tether.perimeter = 24e-3; % [m]
tether.At = (tether.perimeter-tether.thickness*2)/2*tether.thickness; % [m^2]
tether.mt = tether.density*tether.L*tether.At; % tether mass [kg]
[tether.m, tether.phi, tether.LAMBDAt, tether.h_G, tether.Is] = params_2_sbet_params(tether.m1, tether.m2, tether.mt, tether.L); % [kg, radians, non-dimensional, m, kg*m^2]
tether.sigma = 3.5e7; % conductivity for aluminum [1/(Ohm*m)]

%% define constants
constants.mue = 398600; % mu of earth [km^3/s^2]
constants.mue_si = constants.mue*1e9; % mu of earth in meters [m^3/s^2]
A_m2_mum = 8.22e22; % magnetic constant of earths mag field [A*m^2] https://www.britannica.com/science/geomagnetic-field/Dipolar-field
constants.mum = A_m2_mum/795774.71545947673925; % magnetic constant of earths mag field [T*M^3] https://maurermagnetic.com/en/demagnetizing/technology/convert-magnetic-units/
constants.me = 9.1093837e-31; % mass of electron [kg]
constants.mi = 2.6561e-26; % mass of atomic oxygen [kg]
constants.e = 1.60217663e-19; % electron charge [C]
constants.mu = sqrt(constants.me/constants.mi); % Ratio 
constants.gamma = .15e-3; % secondary emission yield [1/V]

%% Define constants
mue = 398600; % mu of earth
mue_si = mue*1e9; % mu of earth in meters
mum = 8e15; % magnetic constant of earths mag field
me = 9.1093837e-31; % mass of electron
mi = 2.6561e-26; % mass of atomic oxygen
e = 1.60217663e-19; % electron charge
mu = sqrt(me/mi); % Ratio (me/mi)^1/2
sigma = 3.5e7; % conductivity for aluminum
gamma = .15e-3; % secondary emission yield

%% load test data
data_02 = load("Em_0_2_1e2_ninf_2e9_2e12_1e2.mat");
nn = length(data_02.Em);
mm = length(data_02.ninf);
test_input = [];
test_target = [];
for ii = 1:nn
    for jj = 1:mm
        if data_02.exitflag(ii,jj) > 0
            test_input(:,end+1) = [tether.L; tether.perimeter; tether.At; data_02.Em(ii); data_02.ninf(jj)];
            test_target(:,end+1) = [data_02.PHIe_0(ii,jj), data_02.h_0(ii,jj)];
        end
    end
end



%% load data
data_04 = load("Em_0_4_1e2_ninf_2e9_2e12_1e2.mat");
data_12 = load("Em_1_2_1e2_ninf_2e9_2e12_1e2.mat");

nn = length(data_04.Em);
mm = length(data_04.ninf);
input = [];
target = [];
for ii = 1:nn
    for jj = 1:mm
        if data_04.exitflag(ii,jj) > 0
            input(:,end+1) = [tether.L; tether.perimeter; tether.At; data_04.Em(ii); data_04.ninf(jj)];
            target(:,end+1) = [data_04.PHIe_0(ii,jj), data_04.h_0(ii,jj)];
        end
    end
end

nn = length(data_12.Em);
mm = length(data_12.ninf);
for ii = 1:nn
    for jj = 1:mm
        if data_12.exitflag(ii,jj) > 0
            input(:,end+1) = [tether.L; tether.perimeter; tether.At; data_12.Em(ii); data_12.ninf(jj)];
            target(:,end+1) = [data_12.PHIe_0(ii,jj), data_12.h_0(ii,jj)];
        end
    end
end

%% normalize
norm_L = tether.L;
norm_perimeter = tether.perimeter;
norm_At = tether.At;
norm_Em = 4;
norm_ninf = 2e12;
norm_phi = 1e5;
norm_h = 1e3;

input(1,:) = input(1,:)/norm_L;
input(2,:) = input(2,:)/norm_perimeter;
input(3,:) = input(3,:)/norm_At;
input(4,:) = input(4,:)/norm_Em;
input(5,:) = input(5,:)/norm_ninf;

target(1,:) = target(1,:)/norm_phi; %phi
target(2,:) = target(2,:)/norm_h; %h

%% normalize test
test_input(1,:) = test_input(1,:)/norm_L;
test_input(2,:) = test_input(2,:)/norm_perimeter;
test_input(3,:) = test_input(3,:)/norm_At;
test_input(4,:) = test_input(4,:)/norm_Em;
test_input(5,:) = test_input(5,:)/norm_ninf;

test_target(1,:) = test_target(1,:)/norm_phi; %phi
test_target(2,:) = test_target(2,:)/norm_h; %h

%% Train
net = feedforwardnet([2].*128,'trainlm');
net = configure(net, input, target);
net.trainParam.show = 50;
net.trainParam.lr = .05;
net.trainParam.epochs = 10000;
net.trainParam.goal = 1e-12;
net.divideParam.trainRatio = .70;
net.divideParam.valRatio = .20;
net.divideParam.testRatio = .10;
net.performFcn = 'mse';
net = train(net, input, target, 'useGPU', 'no');
save("nna_L1N128_1e2_1e2.mat","net")

% load("nna_L4N16_1e2_1e2_100_5km_5mm.mat")

%% show
output = net(input);

% 
error_raw = output - target;
failed_points = 0;
for ii = 1:length(error_raw)
    percent_error_h_0(ii) = 100*abs(error_raw(2,ii))./target(2,ii);
    percent_error_PHIe_0(ii) = 100*abs(error_raw(1,ii))./target(1,ii);
end
for ii = 1:length(error_raw)
    if percent_error_h_0(ii) > 1e4
        percent_error_h_0(ii) = 0;
        failed_points = failed_points + 1;
    end
    if percent_error_PHIe_0(ii) > 1e4
        percent_error_PHIe_0(ii) = 0;
        failed_points = failed_points + 1;
    end
end
fprintf("Performance: %e\n", perform(net, target, net(input)))
fprintf("Percent error training h_0: %f\n", mean(percent_error_h_0))
fprintf("Percent error training PHIe_0: %f\n", mean(percent_error_PHIe_0))

test_output = net(test_input);

test_error_raw = test_output - test_target;
for ii = 1:length(test_error_raw)
    test_percent_error_h_0(ii) = 100*abs(test_error_raw(2,ii))./test_target(2,ii);
    test_percent_error_PHIe_0(ii) = 100*abs(test_error_raw(1,ii))./test_target(1,ii);
end
test_failed_points = 0;
for ii = 1:length(test_error_raw)
    if test_percent_error_h_0(ii) > 1e4
        test_percent_error_h_0(ii) = 0;
        test_failed_points = failed_points + 1;
    end
    if test_percent_error_PHIe_0(ii) > 1e4
        test_percent_error_PHIe_0(ii) = 0;
        test_failed_points = failed_points + 1;
    end
end
fprintf("Performance: %e\n", perform(net, test_target, net(test_input)))
fprintf("Percent error testing h_0: %f\n", mean(test_percent_error_h_0))
fprintf("Percent error testing PHIe_0: %f\n", mean(test_percent_error_PHIe_0))

figure
subplot(2,1,1)
hold on
plot(input(4,:)*norm_Em, target(1,:)*norm_phi,'.')
plot(input(4,:)*norm_Em, output(1,:)*norm_phi,'.')
hold off
xlabel('Em')
ylabel('PHIe_0')
legend('target', 'output')

subplot(2,1,2)
hold on
plot(input(4,:)*norm_Em, target(2,:)*norm_h,'.')
plot(input(4,:)*norm_Em, output(2,:)*norm_h,'.')
hold off
xlabel('Em')
ylabel('h_0')
legend('target', 'output')


figure
subplot(2,1,1)
hold on
plot(input(5,:)*norm_ninf, target(1,:)*norm_phi,'.')
plot(input(5,:)*norm_ninf, output(1,:)*norm_phi,'.')
hold off
xlabel('ninf')
ylabel('PHIe_0')
legend('target', 'output')

subplot(2,1,2)
hold on
plot(input(5,:)*norm_ninf, target(2,:)*norm_h,'.')
plot(input(5,:)*norm_ninf, output(2,:)*norm_h,'.')
hold off
xlabel('ninf')
ylabel('h_0')
legend('target', 'output')


figure
subplot(2,1,1)
hold on
plot(test_input(4,:)*norm_Em, test_target(1,:)*norm_phi,'.')
plot(test_input(4,:)*norm_Em, test_output(1,:)*norm_phi,'.')
hold off
xlabel('Em')
ylabel('PHIe_0')
legend('target', 'output')

subplot(2,1,2)
hold on
plot(test_input(4,:)*norm_Em, test_target(2,:)*norm_h,'.')
plot(test_input(4,:)*norm_Em, test_output(2,:)*norm_h,'.')
hold off
xlabel('Em')
ylabel('h_0')
legend('target', 'output')


figure
subplot(2,1,1)
hold on
plot(test_input(5,:)*norm_ninf, test_target(1,:)*norm_phi,'.')
plot(test_input(5,:)*norm_ninf, test_output(1,:)*norm_phi,'.')
hold off
xlabel('ninf')
ylabel('PHIe_0')
legend('target', 'output')

subplot(2,1,2)
hold on
plot(test_input(5,:)*norm_ninf, test_target(2,:)*norm_h,'.')
plot(test_input(5,:)*norm_ninf, test_output(2,:)*norm_h,'.')
hold off
xlabel('ninf')
ylabel('h_0')
legend('target', 'output')

