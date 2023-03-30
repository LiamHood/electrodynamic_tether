clear; close all; clc;
% save("electron_density.mat");
load("electron_density.mat");
%electrondensity(:,"Jdays"),
tra_input = table2array([electrondensity(:,"LocalTime"),electrondensity(:,"Latitude"), ...
    electrondensity(:,"Longitude"), electrondensity(:,"Altitude")])';
tra_target = table2array(electrondensity(:,"Density"))';

net = feedforwardnet([2,2,2,2,2].^6,'trainlm');
net = configure(net, tra_input, tra_target);
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
% net = train(net, tra_input, tra_target, 'useGPU', 'yes');
net = train(net, tra_input, tra_target, 'useGPU', 'no');
net_tra_output = net(tra_input);

percent_error = 100*abs(net_tra_output-tra_target)./tra_target;
figure
plot(percent_error)

figure
plot(tra_input(1,1:1e3),tra_input(2,1:1e3),'.')
% 
% figure
% plot(tra_input(2,1:1e3))
