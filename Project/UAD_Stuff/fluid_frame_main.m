clc; clear; close all;

%% Definitions
c = 5;
N = 10;

u = 20;
alp = 0;
rho = 1.225;

sn = sind(alp);
cs = cosd(alp);

dx = c/N;
ele_beg = (0:N-1)'*dx;
plate = (0:N)'*dx;

vor_loc = ele_beg + 0.25*dx;
colc_loc = ele_beg + 0.75*dx;

%% Finding Influence Coefficient Matrix
A = zeros(N);
B = zeros(N,1);

% Aij corresponds to ind. vel. at ith collocation point due to jth vortex element
for i = 1:N
    for j = 1:N
        [~,A(i,j)] = getIndVel(1,colc_loc(i),0,vor_loc(j),0);
    end
    [~,B(i)] = getIndVel(1,colc_loc(i),0,c+0.1*dx,0);
end

% Coefficient matrix
C = [A,B;ones(1,N+1)];

%% Time Marching
dt = 0.01;
% dt = 0.01*c/u;
times = (0:dt:5)';

% Flapping Parameters
n = 1;
f1 = 1;
f2 = n*f1;
a = 1;

% Initialisation of important matrices
gbm = zeros(N,length(times));
gw = zeros(size(times));
xw = zeros(size(times));
yw = zeros(size(times));
L = zeros(size(times));
D = zeros(size(times));
vindw = zeros(N,1); 

Bj = zeros(N,length(times)-1);
R = zeros(N+1,1);

status ='o';

% Initialisation of figure for wake capture
fig1 = figure(1);
fig1.Position = [200,200,900,600];
gifname = 'Images\WC_n(8)_f1(1)';
wakeplotname = 'Images\Wake_n(4)_f1(1)';
plot1name = 'Images\Load_n(8)_f1(1)';
plot2name = 'Images\CL_n(8)_f1(1)';

% if isfile([pwd,'\',gifname,'.gif'])
%     fprintf('File Detected!!\n')
%     delete([pwd,'\',gifname,'.gif'])
% end

% For loop for time marching
for m = 1:length(times) % time step index
    disp(m)
    ydot = hdot(times(m), f1, f2, a);

    extflow = u*sn - ydot*cs;
    
    % Origin location at time t
    xo = -u*times(m);
    yo = h(times(m),f1,f2,a);

    % New wake location
    xw(m) = xo + (c + 0.1*dx)*cs;
    yw(m) = yo - (c + 0.1*dx)*sn;
    
    % Current Vortex location
    vor_locx = vor_loc*cs + xo;
    vor_locy = -vor_loc*sn + yo;
    
    % Full plate coordinates
    platex = xo + plate*cs;
    platey = yo - plate*sn;
    
    % Collocation point location
    colc_locx = colc_loc*cs + xo;
    colc_locy = -colc_loc*sn + yo;

    % B vectors and R (RHS) vector calculations
    for i = 1:N
        for j = 1:m-1
            [t1,t2] = getIndVel(1,colc_locx(i),colc_locy(i),xw(j),yw(j));
            Bj(i,j) = t1*sn + t2*cs;
        end
        R(i) = -extflow - sum(Bj(i,1:m-1)*gw(1:m-1)); 
    end
    
    R(end) = -sum(gw(1:m-1));
    
    % Computation of unknown gbm and gwm
    U = C\R;

    gbm(:,m) = U(1:N);
    gw(m) = U(end);

    % Wake Position Update
    for p = 1:m % Position of wake
        [t1, t2] = getIndVel(gbm(:,m), xw(p), yw(p), vor_locx, vor_locy);       % Due to body vortices on wakes
        [t3, t4] = getIndVel(gw(setdiff(1:m,p)), xw(p), yw(p), xw(setdiff(1:m,p)), yw(setdiff(1:m,p)));         % Other wakes on TE wake

        uw = sum(t1) + sum(t3);
        vw = sum(t2) + sum(t4);
        xw(p) = xw(p) + uw*dt;
        yw(p) = yw(p) + vw*dt;
    end

    % Aerodynamic Load Calculations
    if m == 1
        dgbm_dt = sum(gbm(:,m)) / dt;
    else
        dgbm_dt = (sum(gbm(:,m)) - sum(gbm(:,m-1))) / dt;
    end
    
    for i = 1:N
        [t1,t2] = getIndVel(gw(1:m), colc_locx(i), colc_locy(i), xw(1:m), yw(1:m));
        vindw(i) = sum(t1)*sn + sum(t2)*cs;
    end

    L(m) = rho*(u*sum(gbm(:,m)) + dgbm_dt*c);
    D(m) = rho*(vindw'*gbm(:,m) + dgbm_dt*c*deg2rad(alp));


    % Plotting
    
    scatter(xo,yo,10,'red','filled')
    hold on
    plot(platex,platey,'Color','black')
    scatter(xw(1:m),yw(1:m),5,'b','filled')
    hold off
    grid on
    grid minor
    title(['m = ',num2str(m),', t = ',num2str(times(m),'%.3f')])
    xlim([-110,10*c])
    ylim([-10,10])
    if (mod(m,2) == 0) && (status == 'r')
        exportgraphics(gcf,[gifname,'.gif'],'Resolution',150,'Append',true);
    end
    pause(0.01)
end
if status == 'p'
    exportgraphics(gcf,[wakeplotname,'.png'],'Resolution',300)
end


Cl = L/(0.5*rho*u^2*c);
Cd = D/(0.5*rho*u^2*c);

fprintf('Mean CL = %.4f\n',mean(Cl));
fprintf('Mean CD = %.4f\n',mean(Cd));

fig2 = figure(2);
fig2.Position = [200,200,900,600];
plot(times,L,'Color','red','LineWidth',1,'DisplayName','Lift')
hold on
plot(times,D,'Color','blue','LineWidth',1,'DisplayName','Drag')
grid on
grid minor
xlabel('Time t (in s)','FontSize',12)
ylabel('Load (in N)','FontSize',12)
legend('Location','best','FontSize',12)
if status =='r'
    exportgraphics(gcf,[plot1name,'.png'],'Resolution',300);
end

fig3 = figure(3);
fig3.Position = [200,200,900,600];
plot(times,Cl,'Color','red','LineWidth',1,'DisplayName','C_l')
hold on
plot(times,Cd,'Color','blue','LineWidth',1,'DisplayName','C_d')
grid on
grid minor
xlabel('Time t (in s)','FontSize',12)
ylabel('Load Coefficients','FontSize',12)
legend('Location','best','FontSize',12)
if status =='r'
    exportgraphics(gcf,[plot2name,'.png'],'Resolution',300);
end