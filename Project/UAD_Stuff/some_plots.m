clc;
clear;
close all;

t = linspace(0,3,1000);
n = 15;
fs = 0.5;
f1 = 0.5;
f2 = n*f1;
a = 1;
yasm = h(t, f1, f2, a);
ysym = h(t, fs, fs, a);

% f = figure(1);
% f.Position = [200,200,900,600];
% plot(t,ysym)
% hold on
% xline(0)
% yline(0)
% xlim([-0.1,3])
% grid on
% grid('minor')

f = figure(2);
f.Position = [200,200,900,600];
plot(t,yasm)
hold on
xline(0)
yline(0)
xlim([-0.1,0.2])
grid on
grid('minor')