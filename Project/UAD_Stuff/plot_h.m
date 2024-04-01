clc; clear; close all;

% Flapping Parameters
f1 = 1;
a = 1;

nsym = 1;
f2sym = nsym*f1;

nasym = 3;
f2asym = nasym*f1;


T1 = 1/f1;
T2sym = 1/f2sym;
T2asym = 1/f2asym;

T = 0.5*(T1 + T2sym);
times = 0:0.001:T;

ys = h(times,f1,f2sym,a);
yas = h(times,f1,f2asym,a);

fig1 = figure(1);
fig1.Position = [200,200,900,600];
plotname = 'asymvsym';

% if isfile([pwd,'\',plotname,'.png'])
%     fprintf('File Detected!!\n')
%     delete([pwd,'\',plotname,'.png'])
% end

plot(times,ys,'Color','blue','LineWidth',1,'LineStyle','--','DisplayName','Symmetric')
hold on
plot(times,yas,'Color','red','LineStyle','-.','LineWidth',1,'DisplayName','Asymmetric')
grid on
grid minor
xlabel('Time t (in s)','FontSize',12,'FontName','Lucida Fax')
ylabel('h (in m)','FontSize',12,'FontName','Lucida Fax')
xlim([-0.02,1])
xline(0,'Color','black','LineWidth',2.5,'HandleVisibility','off')
yline(0,'Color','black','LineWidth',2.5,'HandleVisibility','off')
title('Symmetric & Asymmetric Flapping','FontSize',12,'FontName','Lucida Fax')
legend('Location','southeast','FontSize',10,'FontName','Lucida Fax')
% exportgraphics(gcf,['asymvsym','.png'],'Resolution',300);


% fig2 = figure(2);
% fig2.Position = [200,200,900,600];
% plot(times,ys,'Color','blue','LineWidth',1,'LineStyle','--')
% grid on
% grid minor
% xlabel('Time t (in s)','FontSize',12,'FontName','Lucida Fax')
% ylabel('h (in m)','FontSize',12,'FontName','Lucida Fax')
% xlim([-0.02,1])
% xline(0,'Color','black','LineWidth',2.5)
% yline(0,'Color','black','LineWidth',2.5)
% title('Symmetric Flapping','FontSize',12,'FontName','Lucida Fax')
% % exportgraphics(gcf,[plot2name,'.png'],'Resolution',300);
% 
% fig3 = figure(3);
% fig3.Position = [200,200,900,600];
% plot(times,yas,'Color','red','LineWidth',1,'LineStyle','--')
% grid on
% grid minor
% xlabel('Time t (in s)','FontSize',12,'FontName','Lucida Fax')
% ylabel('h (in m)','FontSize',12,'FontName','Lucida Fax')
% title('Asymmetric Flapping','FontSize',12,'FontName','Lucida Fax')
% % exportgraphics(gcf,[plot3name,'.png'],'Resolution',300);