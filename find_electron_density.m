function ninf = find_electron_density(jdate, r, e_density_ann)
    utc_time = datetime(jdate,'ConvertFrom','juliandate');
    lla = eci2lla(r'*1e3,[year(utc_time), month(utc_time), day(utc_time),...
        hour(utc_time), minute(utc_time),second(utc_time)]);
    UT = hour(utc_time) + minute(utc_time)/60 + second(utc_time)/3600;
    if lla(2)>0
        LT = UT + lla(2)*(12/180);
    else
        LT = UT + (lla(2)+360)*(12/180);
    end
    
    lat_norm = lla(1)/90;
    long_norm = lla(2)/180;
    alt_norm = (lla(3)/1e3)/2000; % m to km, then normalize to neural net
    LT_norm = LT/24;
    ninf = e_density_ann([lat_norm; long_norm; alt_norm; LT_norm])*2e12; % ionospheric plasma density 

end