clear; close all; clc;


xmesh = linspace(0,pi,10);
solinit = bvpinit(xmesh,@guess);

sol = bvp4c(@bvpfun, @bcfun, solinit);
plot(sol.x,sol.y,'-o')

function dydx = bvpfun(x,y) % equation being solved
dydx = [y(2)
       -y(1)];
end
%-------------------------------------------
function res = bcfun(ya,yb) % boundary conditions
res = [ya(1)
       yb(1)];
end
%-------------------------------------------
function y = guess(x) % guess at solution behavior
y = [sin(x)
     cos(x)];
end
%-------------------------------------------