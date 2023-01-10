clear
utc = [2019 1 4 12 0 0];
lat = 42.283;
long = -71.35;
alt = 1000;
r_ecef = lla2ecef([lat, long, alt]);
% [XYZ,H,D,I,F] = igrfmagm(1000,42.283,-71.35,decyear(2015,7,4),13);
[XYZ,H,D,I,F] = igrfmagm(alt,lat,long,decyear(utc),13);

[X, Y, Z] = ned2ecefv(XYZ(1),XYZ(2),XYZ(3),lat,long);
r_ecef = [-5762640 -1682738 3156028];
v_ecef = [3832 -4024 4837];

[r_eci, xyz_eci] = ecef2eci(utc, r_ecef, [X, Y, Z]);


mag_ned = norm(XYZ)
mag_ecef = norm([X,Y,Z])
mag_eci = norm(xyz_eci)