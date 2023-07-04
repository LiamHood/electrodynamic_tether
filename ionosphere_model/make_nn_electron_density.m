clear; close all; clc;
% save("electron_density.mat");
load("electron_density.mat");
input = table2array([ electrondensity(:,"Latitude"),electrondensity(:,"Longitude"),...
    electrondensity(:,"Altitude"),electrondensity(:,"LT")])';
target = table2array(electrondensity(:,"Density"))';
target = target.*((1e2)^3); % cm^-3 to m^-3

load("testelectron_density.mat");
test_input = table2array([ testelectrondensity(:,"Latitude"),testelectrondensity(:,"Longitude"),...
    testelectrondensity(:,"Altitude"),testelectrondensity(:,"LT")])';
test_target = table2array(testelectrondensity(:,"Density"))';
test_target = test_target.*((1e2)^3); % cm^-3 to m^-3
for ii = 1:length(test_target) % change 0:360 to -180:180
    if test_input(2,ii) > 180
        test_input(2,ii) = test_input(2,ii) - 360;
    end
end
%% normalization

input(1,:) = input(1,:)/90; % latitude, -90 to 90
input(2,:) = input(2,:)/180; % longitude, -180 to 180
input(3,:) = (input(3,:)-100)/2000; % altitude, data made from 100-2000 km
input(4,:) = input(4,:)/24; % time, 0 to 24 hours
target = target/2e12; % density, max is 1.89e12 1/m^3

test_input(1,:) = test_input(1,:)/90; % latitude, -90 to 90
test_input(2,:) = test_input(2,:)/180; % longitude, -180 to 180
test_input(3,:) = test_input(3,:)/2000; % altitude, data made from 100-2000 km
test_input(4,:) = test_input(4,:)/24; % time, 0 to 24 hours
test_target = test_target/2e12; % density, max is 1.89e12 1/m^3


%% Train
net = feedforwardnet([2,2,2,2,2,2].*11,'trainlm');
net = configure(net, input, target);
% for i=1:net.numLayers
%   if strcmp(net.layers{i}.transferFcn,'tansig')
%     net.layers{i}.transferFcn = 'elliotsig';
%   end
% end
net.trainParam.show = 50;
net.trainParam.lr = .05;
net.trainParam.epochs = 10000;
net.trainParam.goal = 1e-12;
net.divideParam.trainRatio = .70;
net.divideParam.valRatio = .20;
net.divideParam.testRatio = .10;
net.performFcn = 'mse';
% net = train(net, tra_input, tra_target, 'useGPU', 'yes');
net = train(net, input, target, 'useGPU', 'no');
save("nna_llat2densityL6N11.mat","net")

% load("nna_llat2densityL5N5.mat","net")
output = net(input);
test_output = net(test_input);

% 
error_raw = output - target;
percent_error = 100*abs(error_raw)./target;
for ii = 1:length(error_raw)
    percent_error(ii) = 100*abs(error_raw(ii))./target(ii);
end
test_error_raw = test_output - test_target;
test_percent_error = 100*abs(test_error_raw)./test_target;
for ii = 1:length(test_error_raw)
    test_percent_error(ii) = 100*abs(test_error_raw(ii))./test_target(ii);
end

fprintf("Percent error training: %f\n", perform(net, target, net(input)))
fprintf("Percent error training: %f\n", mean(percent_error))
fprintf("Percent error test: %f\n", perform(net, test_target, net(test_input)))
fprintf("Percent error test: %f\n", mean(test_percent_error))

altitude = input(3,:)*2000;
latitude = input(1,:)*90;
test_altitude = test_input(3,:)*2000;
test_latitude = test_input(1,:)*90;

figure
hold on
plot(output,".")
plot(target,".")
hold off
ylabel('training density')
legend('net', 'target')

figure
hold on
plot(latitude, output,".")
plot(latitude, target,".")
hold off
ylabel('training density')
xlabel('training latitude')
legend('net', 'target')

figure
hold on
plot(altitude, output,".")
plot(altitude, target,".")
hold off
ylabel('training density')
xlabel('training altitude')
legend('net', 'target')

figure
hold on
plot(test_output,".")
plot(test_target,".")
hold off
ylabel('test density')
legend('net', 'target')

figure
hold on
plot(test_latitude, test_output,".")
plot(test_latitude, test_target,".")
hold off
ylabel('test density')
xlabel('test latitude')
legend('net', 'target')

figure
hold on
plot(test_altitude, test_output,".")
plot(test_altitude, test_target,".")
hold off
ylabel('test density')
xlabel('test altitude')
legend('net', 'target')