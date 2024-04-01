%  Script to find lift of a airfoil over time using Quasi-Steady Approach
clc; clear; close all;

%% Discretization
dt = 0.001;
times = 0:dt:100;

%% Flow & Airfoil Parameters
u = 20;
c = 1;
b = 0.5*c;
xf = 0.25*c;
rho = 1.225;

%% Pitch - Plunge Parameters
a0 = 0.08;
h0 = 1;
% wa = 0.09*u/c;
% wh = 0.09*u/c;
f = 0.025;
wa = 2*pi*f;
wh = 2*pi*f;
ka = wa*c/u;
kh = wh*c/u;
fprintf('Non-D Freq: ka = %.4f\tkh = %.4f\n',ka,kh)

%% Lift & Moment Coefficient Equations
Cl = 2*pi*(alp(times,a0,wa) + hdot(times,h0,wh)/u + alpdot(times,a0,wa)*(0.75*c - xf)/u);
CmLE = - 0.25*Cl - c*pi*alpdot(times,a0,wa)/(8*u);

L = 0.5*rho*u^2*c*Cl;
MLE = 0.5*rho*u^2*c^2*CmLE;

%% Added Mass Corrections
L1 = rho*pi*b^2*(hddot(times,h0,wh) + (0.5*c - xf)*alpddot(times,a0,wa));
L2 = rho*pi*b^2*u*alpdot(times,a0,wa);
L_Cor = L + L1 + L2;
M_Cor = MLE - L1*(0.5*c - xf) - L2*(0.75*c - xf) - rho*pi*b^4*alpddot(times,a0,wa)/8;

%% Plotting
figure(1)
plot(times,L,'Color','red','LineWidth',1,'DisplayName','W/O Cor')
hold on
plot(times,L_Cor,'Color','blue','LineWidth',1,'DisplayName','With Cor')
grid on
grid minor
xlabel('Time t (in sec)','FontSize',14,'FontName','Lucida Fax')
ylabel('Lift L (in N)','FontSize',14,'FontName','Lucida Fax')
title('Quasi-Steady Load Calculation','FontSize',14,'FontName','Lucida Fax')
legend('Location','best','FontName','Lucida Fax','FontSize',10)

figure(2)
plot(times,MLE,'Color','red','LineWidth',1,'DisplayName','W/O Cor')
hold on
plot(times,M_Cor,'Color','blue','LineWidth',1,'DisplayName','With Cor')
grid on
grid minor
xlabel('Time t (in sec)','FontSize',14,'FontName','Lucida Fax')
ylabel('Moment M (in Nm)','FontSize',14,'FontName','Lucida Fax')
title('Quasi-Steady Moment Calculation','FontSize',14,'FontName','Lucida Fax')
legend('Location','best','FontName','Lucida Fax','FontSize',10)