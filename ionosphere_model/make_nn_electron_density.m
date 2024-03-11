clear; close all; clc;
% save("electron_density.mat");
load("electron_density.mat");
load("electron_density_odd_hours.mat");
input = table2array([electrondensity(:,"Latitude"),electrondensity(:,"Longitude"),...
    electrondensity(:,"Altitude"),electrondensity(:,"LT")])';
input_odd = table2array([electrondensityoddhours(:,"Latitude"),electrondensityoddhours(:,"Longitude"),...
    electrondensityoddhours(:,"Altitude"),electrondensityoddhours(:,"LT")])';
input = [input, input_odd];
target = table2array(electrondensity(:,"Density"))';
target_odd = table2array(electrondensityoddhours(:,"Density"))';
target = [target, target_odd];
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
test_target = test_target/2e12; % densx is 1.89e12 1/m^3


%% Train
for p = 0:6
    net = feedforwardnet([2].^p,'trainlm');
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
    name = "nna_llat2densityL1P" + num2str(p) + ".mat";
    save(name,"net")
end
% 
% for p = 0:6
%     net = feedforwardnet([2,2].^p,'trainlm');
%     net = configure(net, input, target);
%     net.trainParam.show = 50;
%     net.trainParam.lr = .05;
%     net.trainParam.epochs = 10000;
%     net.trainParam.goal = 1e-12;
%     net.divideParam.trainRatio = .70;
%     net.divideParam.valRatio = .20;
%     net.divideParam.testRatio = .10;
%     net.performFcn = 'mse';
%     net = train(net, input, target, 'useGPU', 'no');
%     name = "nna_llat2densityL2P" + num2str(p) + ".mat";
%     save(name,"net")
% end

% for p = 0:5
%     net = feedforwardnet([2,2,2].^p,'trainlm');
%     net = configure(net, input, target);
%     net.trainParam.show = 50;
%     net.trainParam.lr = .05;
%     net.trainParam.epochs = 10000;
%     net.trainParam.goal = 1e-12;
%     net.divideParam.trainRatio = .70;
%     net.divideParam.valRatio = .20;
%     net.divideParam.testRatio = .10;
%     net.performFcn = 'mse';
%     net = train(net, input, target, 'useGPU', 'no');
%     name = "nna_llat2densityL3P" + num2str(p) + ".mat";
%     save(name,"net")
% end
% 
% for p = 0:4
%     net = feedforwardnet([2,2,2,2].^p,'trainlm');
%     net = configure(net, input, target);
%     net.trainParam.show = 50;
%     net.trainParam.lr = .05;
%     net.trainParam.epochs = 10000;
%     net.trainParam.goal = 1e-12;
%     net.divideParam.trainRatio = .70;
%     net.divideParam.valRatio = .20;
%     net.divideParam.testRatio = .10;
%     net.performFcn = 'mse';
%     net = train(net, input, target, 'useGPU', 'no');
%     name = "nna_llat2densityL4P" + num2str(p) + ".mat";
%     save(name,"net")
% end
% 
% for p = 0:3
%     net = feedforwardnet([2,2,2,2,2].^p,'trainlm');
%     net = configure(net, input, target);
%     net.trainParam.show = 50;
%     net.trainParam.lr = .0
%     net.trainParam.epochs = 10000;
%     net.trainParam.goal = 1e-12;
%     net.divideParam.trainRatio = .70;
%     net.divideParam.valRatio = .20;
%     net.divideParam.testRatio = .10;
%     net.performFcn = 'mse';
%     net = train(net, input, target, 'useGPU', 'no');
%     name = "nna_llat2densityL5P" + num2str(p) + ".mat";
%     save(name,"net")
% end
% 
% for p = 3:5
%     net = feedforwardnet([2,2,2,2,2,2].^p,'trainlm');
%     net = configure(net, input, target);
%     net.trainParam.show = 50;
%     net.trainParam.lr = .05;
%     net.trainParam.epochs = 10000;
%     net.trainParam.goal = 1e-12;
%     net.divideParam.trainRatio = .70;
%     net.divideParam.valRatio = .20;
%     net.divideParam.testRatio = .10;
%     net.performFcn = 'mse';
%     net = train(net, input, target, 'useGPU', 'no');
%     name = "nna_llat2densityL6P" + num2str(p) + ".mat";
%     save(name,"net")
% end
% 
% p = 7;
% net = feedforwardnet([2].^p,'trainlm');
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
% name = "nna_llat2densityL1P" + num2str(p) + ".mat";
% save(name,"net")
% 
% p = 6;
% net = feedforwardnet([2,2].^p,'trainlm');
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
% name = "nna_llat2densitP" + num2str(p) + ".mat";
% save(name,"net")
% 
