clear
close all
t = 1:1:40*60*60;
for ii = 1:length(t)
    tsec(ii) = t(ii);
    tmin(ii) = 0;
    thour(ii) = 0;
    tday(ii) = 1;
    while tsec(ii) >= 60
        tsec(ii) = tsec(ii) - 60;
        tmin(ii) = tmin(ii) + 1;
    end
    while tmin(ii) >= 60
        tmin(ii) = tmin(ii) - 60;
        thour(ii) = thour(ii) + 1;
    end
    while thour(ii) >= 24
        thour(ii) = thour(ii) - 24;
        tday(ii) = tday(ii) + 1;
    end
end
figure
plot(t,tsec)
figure
plot(t,tmin)
figure
plot(t,thour)
figure
plot(t,tday)