% Script to find the unsteady lift variation over time for a moving airfoil
clear; close all; clc;

%% Discretization
dt = 0.01;
times = 0:dt:24;

%% Flow & Airfoil Parameters
u = 1;
c = 1;
b = 0.5*c;
xf = 0;
e = (xf/c - 0.25);
rho = 1.225;

%% Pitch - Plunge Parameters
a0 = deg2rad(3);
h0 = 0.1;
% wa = 0.09*u/c;
% wh = 0.09*u/c;
f = 0.1;
t0 = 4;
wa = 2*pi*f;
wh = 2*pi*f;
ka = wa*c/u;
kh = wh*c/u;
phih = 0;
phia = acos(2/3);
fprintf('Non-D Freq: ka = %.4f\tkh = %.4f\n\n',ka,kh)

alphavec = zeros(size(times));
alphadvec = zeros(size(times));
hdvec = zeros(size(times));
cnt = 1;
for t = times
    alphavec(cnt) = a(t,t0,a0,wa,phia);
    alphadvec(cnt) = ad(t,t0,a0,wa,phia);
    hdvec(cnt) = hgtd(t,t0,h0,wh,phih);
    cnt = cnt + 1;
end
% a = @(t)(alp(t,a0,wa));
% ad = @(t)(alpdot(t,a0,wa));
% % a = @(t) 0;
% % ad = @(t) 0;
% hgt = @(t)(h(t,h0,wh));
% hgtd = @(t)(hdot(t,h0,wh));

%% Wagner Function Definition
psi1 = 0.165;
psi2 = 0.335;
eps1 = 0.0455;
eps2 = 0.3;
p0 = phi(0,u,b,psi1,psi2,eps1,eps2);
pd0 = phidot(0,u,b,psi1,psi2,eps1,eps2);

% Integrating d1, d2, d3 & d4 over time for given motion
[~,dvec] = ode45(@(t,dvec)dvecdot(t,dvec,@(t)h(t,h0,wh,phih),@(t)alp(t,a0,wa,phia),eps1,eps2,u,b),times,zeros(4,1));
fprintf('DONE WITH DVEC INTEGRATION\n')

%% Computing Unsteady Lift
% L = zeros(size(times));
% d1 = 0; d2 = 0; d3 = 0; d4 = 0;

d1 = dvec(:,1)';
d2 = dvec(:,2)';
d3 = dvec(:,3)';
d4 = dvec(:,4)';

L = zeros(size(times));
for i = 1:length(times)
    t1 = p0*(u*a(times(i),t0,a0,wa,phia) + hgtd(times(i),t0,h0,wh,phih) + (0.75*c - xf)*ad(times(i),t0,a0,wa,phia));
    t2 = phidot(times(i),u,b,psi1,psi2,eps1,eps2)*(hgt(0,t0,h0,wh,phih) + (0.75*c - xf)*a(0,t0,a0,wa,phia)) - pd0*(hgt(times(i),t0,h0,wh,phih) + (0.75*c - xf)*a(times(i),t0,a0,wa,phia));
    t3 = - psi1*(eps1*u/b)^2*d1(i) - psi2*(eps2*u/b)^2*d2(i);
    t4 = (psi1*eps1*u^2/b)*(1 - eps1*(1 - 2*e))*d3(i) + (psi2*eps2*u^2/b)*(1 - eps2*(1 - 2*e))*d4(i);
    L(i) = rho*pi*u*c*(t1 + t2 + t3 + t4);
end

% for t = times
%     d1 = dvec(times==t,1);
%     d2 = dvec(times==t,2);
%     d3 = dvec(times==t,3);
%     d4 = dvec(times==t,4);
% 
%     t1 = p0*(u*a(t) + hgtd(t) + (0.75*c - xf)*ad(t));
%     t2 = phidot(t,u,b,psi1,psi2,eps1,eps2)*(hgt(0) + (0.75*c - xf)*a(0)) - pd0*(hgt(t) + (0.75*c - xf)*a(t));
%     t3 = - psi1*(eps1*u/b)^2*d1 - psi2*(eps2*u/b)^2*d2;
%     t4 = (psi1*eps1*u^2/b)*(1 - eps1*(1 - 2*e))*d3 + (psi2*eps2*u^2/b)*(1 - eps2*(1 - 2*e))*d4;
%     
%     L(times==t) = rho*pi*u*c*(t1 + t2 + t3 + t4);
% %     L(times==t) = rho*pi*u*c*(t1 - t2 + t3 + t4);
%     i = find(times==t);
%     n = length(times);
%     perc = 100*i/n;
%     if abs(rem(perc,5)) < 1.001*100/n
%         fprintf('Completed %3.0f%%\n',perc)
%     end
% 
% %     d1d = h(t,h0,wh) - (eps1*u/b)*d1; d1 = d1 + dt*d1d;
% %     d2d = h(t,h0,wh) - (eps2*u/b)*d2; d2 = d2 + dt*d2d;
% %     d3d = alp(t,a0,wa) - (eps1*u/b)*d3; d3 = d3 + dt*d3d;
% %     d4d = alp(t,a0,wa) - (eps2*u/b)*d4; d4 = d4 + dt*d4d;
% end


%% Plotting the result
figure(1)
plot(times,L,'LineWidth',1,'Color','red')
grid on
grid minor
xlabel('Time t (in sec)','FontSize',14,'FontName','Lucida Fax')
ylabel('Lift L (in N)','FontSize',14,'FontName','Lucida Fax')
title('Unsteady Lift Variation','FontSize',14,'FontName','Lucida Fax')

fprintf('Load at t = 24s = %.4f N\n',L(end))
figure(2)
plot(times,rad2deg(alphavec))

figure(3)
plot(times,hdvec)

figure(4)
plot(times,alphadvec)

%% Functions
function alpha = a(t,t0,a0,wa,phia)
    if t < t0
        alpha = deg2rad(2);
    else
        alpha = alp(t - t0,a0,wa,phia);
    end
end

function alpd = ad(t,t0,a0,wa,phia)
    if t < t0
        alpd = 0;
    else
        alpd = alpdot(t - t0,a0,wa,phia);
    end
end

function y = hgt(t,t0,h0,wh,phih)
    if t < t0
        y = h0;
    else
        y = h(t - t0,h0,wh,phih);
    end
end

function yd = hgtd(t,t0,h0,wh,phih)
    if t < t0
        yd = 0;
    else
        yd = hdot(t - t0,h0,wh,phih);
    end
end