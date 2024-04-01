clc; clear; close all;

% k = 0:0.01:5;
rho = 1.225;
u = 20;
alp0 = deg2rad(5);
f = 1;
c = 1;
b = 0.5*c;
omega = 2*pi*f;
k = omega*b/u;

k0 = besselk(0,1i*k);
k1 = besselk(1,1i*k);
cofk = k1./(k0 + k1);

% Pitching
pL = 1i + (1i + 2/k)*cofk;
Lfreq = pi*rho*u^2*b*alp0*k*pL;
pM = - 1i + 0.25*k + (1i + 2/k)*cofk;
Mfreq = 0.5*pi*rho*u^2*b^2*k*alp0*pM;
fprintf('Pitching\n')
fprintf('k = %.4f\n',k)
fprintf('Lift Amplitude = %.4f %.4fi\n',real(Lfreq),imag(Lfreq))
fprintf('Moment Amplitude = %.4f %.4fi\n',real(Mfreq),imag(Mfreq))

t = 0:0.01:5;
y = alp0*exp(1i*omega*t);
L = Lfreq*exp(1i*omega*t);
M = Mfreq*exp(1i*omega*t);

% Care only about the real part
figure(1)
plot(t,real(L))
hold on
plot(t,real(y))
grid on
grid minor

figure(2)
plot(t,real(L)./abs(L))
hold on
plot(t,real(y)./abs(y))
grid on
grid minor

figure(3)
plot(t,real(M))
hold on
plot(t,real(y))
grid on
grid minor

figure(4)
plot(t,real(M)./abs(M))
hold on
plot(t,real(y)./abs(y))
grid on
grid minor

figure(5)
plot(t,real(L))
hold on
plot(t,real(M))
grid on
grid minor