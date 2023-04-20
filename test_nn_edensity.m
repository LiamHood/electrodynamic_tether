clear; close all; clc;
% save("electron_density.mat");
load("electron_density.mat");
%electrondensity(:,"Jdays"),
tra_input = table2array([electrondensity(:,"Jdays"), electrondensity(:,"LocalTime"),electrondensity(:,"Latitude"), ...
    electrondensity(:,"Longitude"), electrondensity(:,"Altitude")])';
tra_target = table2array(electrondensity(:,"Density"))';
test_input = table2array([testelectrondensity(:,"Jdays"), testelectrondensity(:,"LocalTime"), ...
    testelectrondensity(:,"Latitude"), testelectrondensity(:,"Longitude"), ...
    testelectrondensity(:,"Altitude")])';
test_target = table2array(testelectrondensity(:,"Density"))';

ti_mean = mean(tra_input');
ti_std = std(tra_input');
ta_mean = mean(tra_target');
ta_std = std(tra_target');
for ii = 1:length(ti_mean)
    tra_input(ii, :) = tra_input(ii, :) - ti_mean(ii);
    tra_input(ii, :) = tra_input(ii, :) ./ ti_std(ii);
    test_input(ii, :) = test_input(ii, :) - ti_mean(ii);
    test_input(ii, :) = test_input(ii, :) ./ ti_std(ii);
end
tra_target = tra_target - ta_mean;
tra_target = tra_target./ ta_std;
test_target = test_target - ta_mean;
test_target = test_target./ ta_std;


% net = feedforwardnet([2,2,2,2,2,2].*21,'trainlm');
% net = configure(net, tra_input, tra_target);
% % for i=1:net.numLayers
% %   if strcmp(net.layers{i}.transferFcn,'tansig')
% %     net.layers{i}.transferFcn = 'elliotsig';
% %   end
% % end
% net.trainParam.show = 50;
% net.trainParam.lr = .05;
% net.trainParam.epochs = 10000;
% net.trainParam.goal = 1e-12;
% net.divideParam.trainRatio = .70;
% net.divideParam.valRatio = .20;
% net.divideParam.testRatio = .10;
% net.performFcn = 'mae';
% % net = train(net, tra_input, tra_target, 'useGPU', 'yes');
% net = train(net, tra_input, tra_target, 'useGPU', 'no');

load("net_density_norm_2.mat")
net_test_output = net(test_input);
% save("net_density_norm_2.mat","net")

% error_raw = tra_target - net_tra_output;
% error_actual = (tra_target - net_tra_output)*ta_std + ta_mean;
% % percent_error = 100*abs(net_tra_output-tra_target)./tra_target;
% for ii = 1:length(error_raw)
%     percent_error(ii) = 100*abs(error_raw(ii)/tra_input(ii));
% end
% figure
% plot(percent_error,".")
% % axis([1, length(percent_error), 0, 100])
fprintf("Percent error training: %f\n", perform(net, tra_target, net(tra_input)))
fprintf("Percent error testing: %f\n", perform(net, test_target, net(test_input)))

% figure
% plot(error_raw)
% figure
% plot(error_actual)
% figure
% plot(tra_input(1,1:1e3),tra_input(2,1:1e3),'.')
% 
figure
hold on
plot(net_test_output,".")
plot(test_target,".")
hold off
