clc; clear; close all;

rho = 1.225;
u = 20;
y0 = 0.1;
f = 1;
c = 1;
b = 0.5*c;
omega = 2*pi*f;
k = omega*b/u;

k0 = besselk(0,kg);
k1 = besselk(1,kg);
cofk = k1./(k0 + k1);

j0 = besselj(0,kg);
j1 = besselj(1,kg);
S = (j0 - 1i*j1)*cofk + 1i*j1;