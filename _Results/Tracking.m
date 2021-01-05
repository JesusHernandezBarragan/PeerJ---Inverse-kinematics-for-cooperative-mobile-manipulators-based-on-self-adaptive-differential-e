clear all
close all
clc

s1 = 'CFPSO';
s2 = 'FPA';
s3 = 'SDE';
s4 = 'KABC';

CFPSO_1 = Get_Data('CFPSO_1');
FPA_1 = Get_Data('FPA_1');
IDE_1 = Get_Data('IDE_1');
KABC_1 = Get_Data('KABC_1');

CFPSO_2 = Get_Data('CFPSO_2');
FPA_2 = Get_Data('FPA_2');
IDE_2 = Get_Data('IDE_2');
KABC_2 = Get_Data('KABC_2');

CFPSO_3 = Get_Data('CFPSO_3');
FPA_3 = Get_Data('FPA_3');
IDE_3 = Get_Data('IDE_3');
KABC_3 = Get_Data('KABC_3');

CFPSO_4 = Get_Data('CFPSO_4'); %-----------------
FPA_4 = Get_Data('FPA_4'); %-----------------
IDE_4 = Get_Data('IDE_4'); %-----------------
KABC_4 = Get_Data('KABC_4'); %-----------------


%% Positon error
f = figure;

subplot(2,2,1)
boxplot([sum(CFPSO_1.e)' sum(FPA_1.e)' sum(IDE_1.e)' sum(KABC_1.e)'],'notch','off','labels',{s1,s2,s3,s4})
ylabel('error value (m)');
title('A) Trajectory 1');

subplot(2,2,2)
boxplot([sum(CFPSO_2.e)' sum(FPA_2.e)' sum(IDE_2.e)' sum(KABC_2.e)'],'notch','off','labels',{s1,s2,s3,s4})
ylabel('error value (m)');
title('B) Trajectory 2');

subplot(2,2,3)
boxplot([sum(CFPSO_3.e)' sum(FPA_3.e)' sum(IDE_3.e)' sum(KABC_3.e)'],'notch','off','labels',{s1,s2,s3,s4})
ylabel('error value (m)');
title('C) Trajectory 3');

subplot(2,2,4)
boxplot([sum(CFPSO_4.e)' sum(FPA_4.e)' sum(IDE_4.e)' sum(KABC_4.e)'],'notch','off','labels',{s1,s2,s3,s4})
ylabel('error value (m)');
title('D) Trajectory 4');

% Save_BoxPlot(f,'PositionError')
Save_BoxPlot(f,'Figure5')

%% Motion error
f = figure;

subplot(2,2,1)
boxplot([CFPSO_1.eq' FPA_1.eq' IDE_1.eq' KABC_1.eq'],'notch','off','labels',{s1,s2,s3,s4})
ylabel('error value');
title('A) Trajectory 1');

subplot(2,2,2)
boxplot([CFPSO_2.eq' FPA_2.eq' IDE_2.eq' KABC_2.eq'],'notch','off','labels',{s1,s2,s3,s4})
ylabel('error value');
title('B) Trajectory 2');

subplot(2,2,3)
boxplot([CFPSO_3.eq' FPA_3.eq' IDE_3.eq' KABC_3.eq'],'notch','off','labels',{s1,s2,s3,s4})
ylabel('error value');
title('C) Trajectory 3');

subplot(2,2,4)
boxplot([CFPSO_4.eq' FPA_4.eq' IDE_4.eq' KABC_4.eq'],'notch','off','labels',{s1,s2,s3,s4})
ylabel('error value');
title('D) Trajectory 4');

% Save_Figure(f,'MotionError')
Save_BoxPlot(f,'Figure6')

%% Execution Time
f = figure;

subplot(2,2,1)
boxplot([CFPSO_1.time' FPA_1.time' IDE_1.time' KABC_1.time'],'notch','off','labels',{s1,s2,s3,s4})
ylabel('time (s)');
title('A) Trajectory 1');

subplot(2,2,2)
boxplot([CFPSO_2.time' FPA_2.time' IDE_2.time' KABC_2.time'],'notch','off','labels',{s1,s2,s3,s4})
ylabel('time (s)');
title('B) Trajectory 2');

subplot(2,2,3)
boxplot([CFPSO_3.time' FPA_3.time' IDE_3.time' KABC_3.time'],'notch','off','labels',{s1,s2,s3,s4})
ylabel('time (s)');
title('C) Trajectory 3');

subplot(2,2,4)
boxplot([CFPSO_4.time' FPA_4.time' IDE_4.time' KABC_4.time'],'notch','off','labels',{s1,s2,s3,s4})
ylabel('time (s)');
title('D) Trajectory 4');

% Save_Figure(f,'ExecutionTime')
Save_BoxPlot(f,'Figure7')


%% Tracking results
f = figure;

subplot(4,2,1)
hold on
grid on
plot3(IDE_1.td1(1,:),IDE_1.td1(2,:),IDE_1.td1(3,:),'-','LineWidth',1.5)
plot3(IDE_1.t1(1,:),IDE_1.t1(2,:),IDE_1.t1(3,:),'--','LineWidth',1.5)
xlim([0.49 0.51])
ylim([0.40 1.60])
zlim([0.34 0.46])
view([-60 20])
xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')
title('A) Trajectory 1')

subplot(4,2,2)
hold on
grid on
plot3(IDE_2.td1(1,:),IDE_2.td1(2,:),IDE_2.td1(3,:),'-','LineWidth',1.5)
plot3(IDE_2.t1(1,:),IDE_2.t1(2,:),IDE_2.t1(3,:),'--','LineWidth',1.5)
xlim([0.49 0.51])
ylim([-0.06 0.06])
zlim([0.34 0.46])
view([-60 20])
xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')
title('B) Trajectory 2')
% legend({'$\mathbf{t}_{e_{1}}^{*}$','$\mathbf{t}_{e_{1}}^{1}$'},'Interpreter','latex','Location','best','Location','best')

subplot(4,2,3)
hold on
grid on
plot3(IDE_3.td1(1,:),IDE_3.td1(2,:),IDE_3.td1(3,:),'-','LineWidth',1.5)
plot3(IDE_3.t1(1,:),IDE_3.t1(2,:),IDE_3.t1(3,:),'--','LineWidth',1.5)
xlim([0.49 0.51])
ylim([0.40 1.60])
zlim([0.34 0.46])
view([-60 20])
xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')
title('C) Trajectory 3')

subplot(4,2,4)
hold on
grid on
plot3(IDE_4.td1(1,:),IDE_4.td1(2,:),IDE_4.td1(3,:),'-','LineWidth',1.5)
plot3(IDE_4.t1(1,:),IDE_4.t1(2,:),IDE_4.t1(3,:),'--','LineWidth',1.5)
xlim([0.49 0.51])
ylim([-0.06 0.06])
zlim([0.34 0.46])
view([-60 20])
xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')
title('D) Trajectory 4')

subplot(4,2,5)
hold on
grid on
plot3(IDE_1.td2(1,:),IDE_1.td2(2,:),IDE_1.td2(3,:),'-','LineWidth',1.5)
plot3(IDE_1.t2(1,:),IDE_1.t2(2,:),IDE_1.t2(3,:),'--','LineWidth',1.5)
xlim([0.49 0.51])
ylim([0.40-1 1.60-1])
zlim([0.34 0.46])
view([-60 20])
xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')
title('E) Trajectory 1')

subplot(4,2,6)
hold on
grid on
plot3(IDE_2.td2(1,:),IDE_2.td2(2,:),IDE_2.td2(3,:),'-','LineWidth',1.5)
plot3(IDE_2.t2(1,:),IDE_2.t2(2,:),IDE_2.t2(3,:),'--','LineWidth',1.5)
xlim([0.49 0.51])
ylim([-0.06-1 0.06-1])
zlim([0.34 0.46])
view([-60 20])
xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')
title('F) Trajectory 2')

subplot(4,2,7)
hold on
grid on
plot3(IDE_3.td2(1,:),IDE_3.td2(2,:),IDE_3.td2(3,:),'-','LineWidth',1.5)
plot3(IDE_3.t2(1,:),IDE_3.t2(2,:),IDE_3.t2(3,:),'--','LineWidth',1.5)
xlim([0.49 0.51])
ylim([0.40-1 1.60-1])
zlim([0.34 0.46])
view([-60 20])
xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')
title('G) Trajectory 3')

subplot(4,2,8)
hold on
grid on
plot3(IDE_4.td2(1,:),IDE_4.td2(2,:),IDE_4.td2(3,:),'-','LineWidth',1.5)
plot3(IDE_4.t2(1,:),IDE_4.t2(2,:),IDE_4.t2(3,:),'--','LineWidth',1.5)
xlim([0.49 0.51])
ylim([-0.06-1 0.06-1])
zlim([0.34 0.46])
view([-60 20])
xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')
title('H) Trajectory 4')

Save_Figure(f,'Figure8',{'$\mathbf{t}^{*}$','$\mathbf{t}$'})


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

f = figure;

subplot(4,2,1)
hold on
grid on
plot3(KABC_1.td1(1,:),KABC_1.td1(2,:),KABC_1.td1(3,:),'-','LineWidth',1.5)
plot3(KABC_1.t1(1,:),KABC_1.t1(2,:),KABC_1.t1(3,:),'--','LineWidth',1.5)
xlim([0.49 0.51])
ylim([0.40 1.60])
zlim([0.34 0.46])
view([-60 20])
xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')
title('A) Trajectory 1')

subplot(4,2,2)
hold on
grid on
plot3(KABC_2.td1(1,:),KABC_2.td1(2,:),KABC_2.td1(3,:),'-','LineWidth',1.5)
plot3(KABC_2.t1(1,:),KABC_2.t1(2,:),KABC_2.t1(3,:),'--','LineWidth',1.5)
xlim([0.49 0.51])
ylim([-0.06 0.06])
zlim([0.34 0.46])
view([-60 20])
xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')
title('B) Trajectory 2')
% legend({'$\mathbf{t}_{e_{1}}^{*}$','$\mathbf{t}_{e_{1}}^{1}$'},'Interpreter','latex','Location','best')

subplot(4,2,3)
hold on
grid on
plot3(KABC_3.td1(1,:),KABC_3.td1(2,:),KABC_3.td1(3,:),'-','LineWidth',1.5)
plot3(KABC_3.t1(1,:),KABC_3.t1(2,:),KABC_3.t1(3,:),'--','LineWidth',1.5)
xlim([0.49 0.51])
ylim([0.40 1.60])
zlim([0.34 0.46])
view([-60 20])
xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')
title('C) Trajectory 3')

subplot(4,2,4)
hold on
grid on
plot3(KABC_4.td1(1,:),KABC_4.td1(2,:),KABC_4.td1(3,:),'-','LineWidth',1.5)
plot3(KABC_4.t1(1,:),KABC_4.t1(2,:),KABC_4.t1(3,:),'--','LineWidth',1.5)
xlim([0.49 0.51])
ylim([-0.06 0.06])
zlim([0.34 0.46])
view([-60 20])
xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')
title('D) Trajectory 4')

subplot(4,2,5)
hold on
grid on
plot3(KABC_1.td2(1,:),KABC_1.td2(2,:),KABC_1.td2(3,:),'-','LineWidth',1.5)
plot3(KABC_1.t2(1,:),KABC_1.t2(2,:),KABC_1.t2(3,:),'--','LineWidth',1.5)
xlim([0.49 0.51])
ylim([0.40-1 1.60-1])
zlim([0.34 0.46])
view([-60 20])
xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')
title('E) Trajectory 1')

subplot(4,2,6)
hold on
grid on
plot3(KABC_2.td2(1,:),KABC_2.td2(2,:),KABC_2.td2(3,:),'-','LineWidth',1.5)
plot3(KABC_2.t2(1,:),KABC_2.t2(2,:),KABC_2.t2(3,:),'--','LineWidth',1.5)
xlim([0.49 0.51])
ylim([-0.06-1 0.06-1])
zlim([0.34 0.46])
view([-60 20])
xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')
title('F) Trajectory 2')

subplot(4,2,7)
hold on
grid on
plot3(KABC_3.td2(1,:),KABC_3.td2(2,:),KABC_3.td2(3,:),'-','LineWidth',1.5)
plot3(KABC_3.t2(1,:),KABC_3.t2(2,:),KABC_3.t2(3,:),'--','LineWidth',1.5)
xlim([0.49 0.51])
ylim([0.40-1 1.60-1])
zlim([0.34 0.46])
view([-60 20])
xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')
title('G) Trajectory 3')

subplot(4,2,8)
hold on
grid on
plot3(KABC_4.td2(1,:),KABC_4.td2(2,:),KABC_4.td2(3,:),'-','LineWidth',1.5)
plot3(KABC_4.t2(1,:),KABC_4.t2(2,:),KABC_4.t2(3,:),'--','LineWidth',1.5)
xlim([0.49 0.51])
ylim([-0.06-1 0.06-1])
zlim([0.34 0.46])
view([-60 20])
xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')
title('H) Trajectory 4')

Save_Figure(f,'Figure9',{'$\mathbf{t}^{*}$','$\mathbf{t}$'})



%% Funciones 
function data = Get_Data (s)
    load(s)

    Txy = @(d) [1 0 0 d; 0 1 0 0; 0 0 1 0; 0 0 0 1];
    Tyb = @(d) [1 0 0 0; 0 1 0 d; 0 0 1 0; 0 0 0 1];
    Tb0 = @(theta) [cos(theta), -sin(theta) 0, 0.140; sin(theta), cos(theta), 0, 0; ...
                    0, 0, 1.0, 0.151; 0, 0, 0, 1.0];
    T01 = @(theta) [cos(theta), 0, sin(theta), 0.033*cos(theta); sin(theta), 0, -cos(theta), 0.033*sin(theta); ...
                    0, 1.0, 0, 0.147; 0, 0, 0, 1.0];
    T12 = @(theta) [cos(theta), -sin(theta), 0, 0.155*cos(theta); sin(theta), cos(theta), 0, 0.155*sin(theta); ...
                    0, 0, 1.0, 0; 0, 0, 0, 1.0];
    T23 = @(theta) [cos(theta), -sin(theta), 0, 0.135*cos(theta); sin(theta), cos(theta), 0, 0.135*sin(theta); ...
                    0, 0, 1.0, 0; 0, 0, 0, 1.0];
    T34 = @(theta) [cos(theta), 0, sin(theta), 0; sin(theta), 0, -cos(theta), 0; ...
                    0, 1.0, 0, 0; 0, 0, 0, 1.0];
    T45 = @(theta) [cos(theta), -sin(theta), 0, 0; sin(theta), cos(theta), 0, 0; ...
                    0, 0, 1.0, 0.2174; 0, 0, 0, 1.0];

    [~,N] = size(Q);

    e = zeros(2,N);
    eq = zeros(1,N);
    td1 = zeros(3,N);
    td2 = zeros(3,N);
    t1 = zeros(3,N);
    t2 = zeros(3,N);

    q0 = [0.0 0.5 pi/4 0 pi/2 -pi/4 pi/4 0.0 0.0 -0.5 pi/4 0 pi/2 -pi/4 pi/4 0.0]'; 
    
    for i=1:N
        td1(:,i) = [x(i) y(i) z(i)]';
        td2(:,i) = td1(:,i) + prd;

        xgb = Q(:,i);

        q1 = xgb(1:8);
        Tx = Txy(q1(1));
        Ty = Tx*Tyb(q1(2));
        Tb = Ty*Tb0(q1(3));
        T1 = Tb*T01(q1(4));
        T2 = T1*T12(q1(5));
        T3 = T2*T23(q1(6));
        T4 = T3*T34(q1(7));
        T5 = T4*T45(q1(8));
        Tq1 = T5;

        q2 = xgb(9:16);
        Tx = Txy(q2(1));
        Ty = Tx*Tyb(q2(2));
        Tb = Ty*Tb0(q2(3));
        T1 = Tb*T01(q2(4));
        T2 = T1*T12(q2(5));
        T3 = T2*T23(q2(6));
        T4 = T3*T34(q2(7));
        T5 = T4*T45(q2(8));
        Tq2 = T5;

        t1(:,i) = Tq1(1:3,4);
        t2(:,i) = Tq2(1:3,4);

        e(1,i) = norm(td1(:,i)-t1(:,i));
        e(2,i) = norm(td2(:,i)-t2(:,i));
        
        if i~=N
            eq(i+1) = norm(Q(:,i+1)-Q(:,i));
        else
            eq(1) = norm(q0-Q(:,1));
        end
    end

    data.e = e;
    data.eq = eq;
    data.time = time;
    data.N = N;
    data.fit = fitness;
    data.Q = Q;
    data.td1 = td1;
    data.td2 = td2;
    data.t1 = t1;
    data.t2 = t2;
end

function Save_BoxPlot (f,name)
    ha = get(gcf,'children');
    
    set(ha(4),'position',[0.09 0.57 0.39 0.37]);
    set(ha(2),'position',[0.09 0.07 0.39 0.37]);
    
    set(ha(3),'position',[0.59 0.57 0.39 0.37]);
    set(ha(1),'position',[0.59 0.07 0.39 0.37]);

    print(f,['TrackingResults/' name '.png'],'-dpng','-r300');
%     print(f,['TrackingResults/' name '.eps'],'-depsc','-r300');
end

function Save_Figure (f,name,slegend)
    ax = gcf;
    set (gcf, 'Position', [ax.Position(1)*0.1 ax.Position(2)*0.1 ax.Position(3) ax.Position(4)*2])

    ha = get(gcf,'children');
    
    set(ha(2),'position',[0.00+0.12 0.00+0.05 0.30 0.17]);
    set(ha(1),'position',[0.50+0.12 0.00+0.05 0.30 0.17]);
    
    set(ha(4),'position',[0.00+0.12 0.25+0.05 0.30 0.17]);
    set(ha(3),'position',[0.50+0.12 0.25+0.05 0.30 0.17]);

    set(ha(6),'position',[0.00+0.12 0.50+0.05 0.30 0.17]);
    set(ha(5),'position',[0.50+0.12 0.50+0.05 0.30 0.17]);

    set(ha(8),'position',[0.00+0.12 0.75+0.05 0.30 0.17]);
    set(ha(7),'position',[0.50+0.12 0.75+0.05 0.30 0.17]);

    legend(ha(7),slegend,'Location','best','Interpreter','latex')
    
    print(f,['TrackingResults/' name '.png'],'-dpng','-r300');
%     print(f,['TrackingResults/' name '.eps'],'-depsc','-r300');
end

%     set(ha(1),'position',[.5 .1 .4 .4])
%     set(ha(2),'position',[.1 .1 .4 .4])
%     set(ha(3),'position',[.5 .5 .4 .4])
%     set(ha(4),'position',[.1 .5 .4 .4])

%     ax = gca; outerpos = ax.OuterPosition; ti = ax.TightInset;  
%     left = outerpos(1) + ti(1); bottom = outerpos(2) + ti(2); 
%     ax_width = outerpos(3) - ti(1) - ti(3); ax_height = outerpos(4) - ti(2) - ti(4); 
%     ax.Position = [left bottom ax_width ax_height];

% 	set (gcf, 'Position', [300 300 700 500])
