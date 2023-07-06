clear; close all; clc;
% https://arc-aiaa-org.proxy.wichita.edu/doi/pdf/10.2514/6.2004-5309

%% Define tether parameters
L = 5000; % tether length
tether_radius = .0005;
per = tether_radius*2*pi; %perimeter in meters
At = tether_radius^2*pi;

%% load test data
data_10 = load("starting_values_1e2_1e2_10_5km_5mm.mat");
nn = length(data_10.Em);
mm = length(data_10.ninf);
test_input = [];
test_target = [];
for ii = 1:nn
    for jj = 1:mm
        if data_10.exitflag(ii,jj) == 1
            test_input(:,end+1) = [L; per; At; data_10.Em(ii); data_10.ninf(jj)];
            test_target(:,end+1) = [data_10.h_0(ii,jj), data_10.I_B(ii,jj), data_10.PHIe_0(ii,jj)];
        end
    end
end



%% load data
data_100 = load("starting_values_1e2_1e2_100_5km_5mm.mat");
data_1 = load("starting_values_1e2_1e2_1_5km_5mm.mat");

nn = length(data_100.Em);
mm = length(data_100.ninf);
input = [];
target = [];
for ii = 1:nn
    for jj = 1:mm
        if data_100.exitflag(ii,jj) == 1
            input(:,end+1) = [L; per; At; data_100.Em(ii); data_100.ninf(jj)];
            target(:,end+1) = [data_100.h_0(ii,jj), data_100.I_B(ii,jj), data_100.PHIe_0(ii,jj)];
        end
    end
end

nn = length(data_1.Em);
mm = length(data_1.ninf);
for ii = 1:nn
    for jj = 1:mm
        if data_1.exitflag(ii,jj) == 1
            input(:,end+1) = [L; per; At; data_1.Em(ii); data_1.ninf(jj)];
            target(:,end+1) = [data_1.h_0(ii,jj), data_1.I_B(ii,jj), data_1.PHIe_0(ii,jj)];
        end
    end
end

%% normalize
input(1,:) = input(1,:)/L;
input(2,:) = input(2,:)/per;
input(3,:) = input(3,:)/At;
input(4,:) = input(4,:)/max(data_100.Em);
input(5,:) = input(5,:)/max(data_100.ninf);

target(1,:) = target(1,:)/1e3;
target(2,:) = target(2,:)/50;
target(3,:) = target(3,:)/1e5;

%% normalize test
test_input(1,:) = test_input(1,:)/L;
test_input(2,:) = test_input(2,:)/per;
test_input(3,:) = test_input(3,:)/At;
test_input(4,:) = test_input(4,:)/max(data_100.Em);
test_input(5,:) = test_input(5,:)/max(data_100.ninf);

test_target(1,:) = test_target(1,:)/1e3;
test_target(2,:) = test_target(2,:)/50;
test_target(3,:) = test_target(3,:)/1e5;

%% Train
% net = feedforwardnet([2,2,2,2].*32,'trainlm');
% net = configure(net, input, target);
% net.trainParam.show = 50;
% net.trainParam.lr = .05;
% net.trainParam.epochs = 10000;
% net.trainParam.goal = 1e-12;
% net.divideParam.trainRatio = .70;
% net.divideParam.valRatio = .20;
% net.divideParam.testRatio = .10;
% net.performFcn = 'mse';
% net = train(net, input, target, 'useGPU', 'no');
% save("nna_L4N16_1e2_1e2_100_5km_5mm.mat","net")
load("nna_L4N16_1e2_1e2_100_5km_5mm.mat")
%% show
output = net(input);

% 
error_raw = output - target;
failed_points = 0;
for ii = 1:length(error_raw)
    percent_error_h_0(ii) = 100*abs(error_raw(1,ii))./target(1,ii);
    percent_error_I_B(ii) = 100*abs(error_raw(2,ii))./target(2,ii);
    percent_error_PHIe_0(ii) = 100*abs(error_raw(3,ii))./target(3,ii);
end
for ii = 1:length(error_raw)
    if percent_error_h_0(ii) > 1e4
        percent_error_h_0(ii) = 0;
        failed_points = failed_points + 1;
    end
    if percent_error_I_B(ii) > 1e4
        percent_error_I_B(ii) = 0;
        failed_points = failed_points + 1;
    end
    if percent_error_PHIe_0(ii) > 1e4
        percent_error_PHIe_0(ii) = 0;
        failed_points = failed_points + 1;
    end
end
fprintf("Performance: %e\n", perform(net, target, net(input)))
fprintf("Percent error training h_0: %f\n", mean(percent_error_h_0))
fprintf("Percent error training I_B: %f\n", mean(percent_error_I_B))
fprintf("Percent error training PHIe_0: %f\n", mean(percent_error_PHIe_0))

test_output = net(test_input);

test_error_raw = test_output - test_target;
for ii = 1:length(test_error_raw)
    test_percent_error_h_0(ii) = 100*abs(test_error_raw(1,ii))./test_target(1,ii);
    test_percent_error_I_B(ii) = 100*abs(test_error_raw(2,ii))./test_target(2,ii);
    test_percent_error_PHIe_0(ii) = 100*abs(test_error_raw(3,ii))./test_target(3,ii);
end
test_failed_points = 0;
for ii = 1:length(test_error_raw)
    if test_percent_error_h_0(ii) > 1e4
        test_percent_error_h_0(ii) = 0;
        test_failed_points = failed_points + 1;
    end
    if test_percent_error_I_B(ii) > 1e4
        test_percent_error_I_B(ii) = 0;
        test_failed_points = failed_points + 1;
    end
    if test_percent_error_PHIe_0(ii) > 1e4
        test_percent_error_PHIe_0(ii) = 0;
        test_failed_points = failed_points + 1;
    end
end
fprintf("Performance: %e\n", perform(net, test_target, net(test_input)))
fprintf("Percent error testing h_0: %f\n", mean(test_percent_error_h_0))
fprintf("Percent error testing I_B: %f\n", mean(test_percent_error_I_B))
fprintf("Percent error testing PHIe_0: %f\n", mean(test_percent_error_PHIe_0))

figure
subplot(3,1,1)
hold on
plot(input(4,:)*100, target(1,:)*1e3,'.')
plot(input(4,:)*100, output(1,:)*1e3,'.')
hold off
xlabel('Em')
ylabel('h_0')
legend('target', 'output')

subplot(3,1,2)
hold on
plot(input(4,:)*100, target(2,:)*50,'.')
plot(input(4,:)*100, output(2,:)*50,'.')
hold off
xlabel('Em')
ylabel('I_B')
legend('target', 'output')

subplot(3,1,3)
hold on
plot(input(4,:)*100, target(3,:)*1e5,'.')
plot(input(4,:)*100, output(3,:)*1e5,'.')
hold off
xlabel('Em')
ylabel('PHIe_0')
legend('target', 'output')

figure
subplot(3,1,1)
hold on
plot(input(5,:)*max(data_100.ninf), target(1,:)*1e3,'.')
plot(input(5,:)*max(data_100.ninf), output(1,:)*1e3,'.')
hold off
xlabel('ninf')
ylabel('h_0')
legend('target', 'output')

subplot(3,1,2)
hold on
plot(input(5,:)*max(data_100.ninf), target(2,:)*50,'.')
plot(input(5,:)*max(data_100.ninf), output(2,:)*50,'.')
hold off
xlabel('ninf')
ylabel('I_B')
legend('target', 'output')

subplot(3,1,3)
hold on
plot(input(5,:)*max(data_100.ninf), target(3,:)*1e5,'.')
plot(input(5,:)*max(data_100.ninf), output(3,:)*1e5,'.')
hold off
xlabel('ninf')
ylabel('PHIe_0')
legend('target', 'output')

figure
subplot(3,1,1)
hold on
plot(test_input(4,:)*100, test_target(1,:)*1e3,'.')
plot(test_input(4,:)*100, test_output(1,:)*1e3,'.')
hold off
xlabel('Em')
ylabel('h_0')
legend('target', 'output')

subplot(3,1,2)
hold on
plot(test_input(4,:)*100, test_target(2,:)*50,'.')
plot(test_input(4,:)*100, test_output(2,:)*50,'.')
hold off
xlabel('Em')
ylabel('I_B')
legend('target', 'output')

subplot(3,1,3)
hold on
plot(test_input(4,:)*100, test_target(3,:)*1e5,'.')
plot(test_input(4,:)*100, test_output(3,:)*1e5,'.')
hold off
xlabel('Em')
ylabel('PHIe_0')
legend('target', 'output')

figure
subplot(3,1,1)
hold on
plot(test_input(5,:)*max(data_100.ninf), test_target(1,:)*1e3,'.')
plot(test_input(5,:)*max(data_100.ninf), test_output(1,:)*1e3,'.')
hold off
xlabel('ninf')
ylabel('h_0')
legend('target', 'output')

subplot(3,1,2)
hold on
plot(test_input(5,:)*max(data_100.ninf), test_target(2,:)*50,'.')
plot(test_input(5,:)*max(data_100.ninf), test_output(2,:)*50,'.')
hold off
xlabel('ninf')
ylabel('I_B')
legend('target', 'output')

subplot(3,1,3)
hold on
plot(test_input(5,:)*max(data_100.ninf), test_target(3,:)*1e5,'.')
plot(test_input(5,:)*max(data_100.ninf), test_output(3,:)*1e5,'.')
hold off
xlabel('ninf')
ylabel('PHIe_0')
legend('target', 'output')