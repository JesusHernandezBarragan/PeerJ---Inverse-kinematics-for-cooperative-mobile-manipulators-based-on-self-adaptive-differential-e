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

%% Joint Motion 
f = figure;

subplot(4,2,1)
hold on
grid on
plot(IDE_1.Q(1,:),'LineWidth',1.5)
plot(IDE_1.Q(2,:),'LineWidth',1.5)
plot(IDE_1.Q(3,:),'LineWidth',1.5)
plot(IDE_1.Q(4,:),'LineWidth',1.5)
plot(IDE_1.Q(5,:),'LineWidth',1.5)
plot(IDE_1.Q(6,:),'LineWidth',1.5)
plot(IDE_1.Q(7,:),'LineWidth',1.5)
plot(IDE_1.Q(8,:),'LineWidth',1.5,'Color',[0 0 1])
xlabel('k-th trajectory point')
ylabel('joint value (m,rad)')
title('A) Trajectory 1')

subplot(4,2,2)
hold on
grid on
plot(IDE_2.Q(1,:),'LineWidth',1.5)
plot(IDE_2.Q(2,:),'LineWidth',1.5)
plot(IDE_2.Q(3,:),'LineWidth',1.5)
plot(IDE_2.Q(4,:),'LineWidth',1.5)
plot(IDE_2.Q(5,:),'LineWidth',1.5)
plot(IDE_2.Q(6,:),'LineWidth',1.5)
plot(IDE_2.Q(7,:),'LineWidth',1.5)
plot(IDE_2.Q(8,:),'LineWidth',1.5,'Color',[0 0 1])
xlabel('k-th trajectory point')
ylabel('joint value (m,rad)')
title('B) Trajectory 2')

subplot(4,2,3)
hold on
grid on
plot(IDE_3.Q(1,:),'LineWidth',1.5)
plot(IDE_3.Q(2,:),'LineWidth',1.5)
plot(IDE_3.Q(3,:),'LineWidth',1.5)
plot(IDE_3.Q(4,:),'LineWidth',1.5)
plot(IDE_3.Q(5,:),'LineWidth',1.5)
plot(IDE_3.Q(6,:),'LineWidth',1.5)
plot(IDE_3.Q(7,:),'LineWidth',1.5)
plot(IDE_3.Q(8,:),'LineWidth',1.5,'Color',[0 0 1])
xlabel('k-th trajectory point')
ylabel('joint value (m,rad)')
title('C) Trajectory 3')

subplot(4,2,4)
hold on
grid on
plot(IDE_4.Q(1,:),'LineWidth',1.5)
plot(IDE_4.Q(2,:),'LineWidth',1.5)
plot(IDE_4.Q(3,:),'LineWidth',1.5)
plot(IDE_4.Q(4,:),'LineWidth',1.5)
plot(IDE_4.Q(5,:),'LineWidth',1.5)
plot(IDE_4.Q(6,:),'LineWidth',1.5)
plot(IDE_4.Q(7,:),'LineWidth',1.5)
plot(IDE_4.Q(8,:),'LineWidth',1.5,'Color',[0 0 1])
xlabel('k-th trajectory point')
ylabel('joint value (m,rad)')
title('D) Trajectory 4')

subplot(4,2,5)
hold on
grid on
plot(IDE_1.Q(9,:),'LineWidth',1.5)
plot(IDE_1.Q(10,:),'LineWidth',1.5)
plot(IDE_1.Q(11,:),'LineWidth',1.5)
plot(IDE_1.Q(12,:),'LineWidth',1.5)
plot(IDE_1.Q(13,:),'LineWidth',1.5)
plot(IDE_1.Q(14,:),'LineWidth',1.5)
plot(IDE_1.Q(15,:),'LineWidth',1.5)
plot(IDE_1.Q(16,:),'LineWidth',1.5,'Color',[0 0 1])
xlabel('k-th trajectory point')
ylabel('joint value (m,rad)')
title('E) Trajectory 1')

subplot(4,2,6)
hold on
grid on
plot(IDE_2.Q(9,:),'LineWidth',1.5)
plot(IDE_2.Q(10,:),'LineWidth',1.5)
plot(IDE_2.Q(11,:),'LineWidth',1.5)
plot(IDE_2.Q(12,:),'LineWidth',1.5)
plot(IDE_2.Q(13,:),'LineWidth',1.5)
plot(IDE_2.Q(14,:),'LineWidth',1.5)
plot(IDE_2.Q(15,:),'LineWidth',1.5)
plot(IDE_2.Q(16,:),'LineWidth',1.5,'Color',[0 0 1])
xlabel('k-th trajectory point')
ylabel('joint value (m,rad)')
title('F) Trajectory 2')

subplot(4,2,7)
hold on
grid on
plot(IDE_3.Q(9,:),'LineWidth',1.5)
plot(IDE_3.Q(10,:),'LineWidth',1.5)
plot(IDE_3.Q(11,:),'LineWidth',1.5)
plot(IDE_3.Q(12,:),'LineWidth',1.5)
plot(IDE_3.Q(13,:),'LineWidth',1.5)
plot(IDE_3.Q(14,:),'LineWidth',1.5)
plot(IDE_3.Q(15,:),'LineWidth',1.5)
plot(IDE_3.Q(16,:),'LineWidth',1.5,'Color',[0 0 1])
xlabel('k-th trajectory point')
ylabel('joint value (m,rad)')
title('G) Trajectory 3')

subplot(4,2,8)
hold on
grid on
plot(IDE_4.Q(9,:),'LineWidth',1.5)
plot(IDE_4.Q(10,:),'LineWidth',1.5)
plot(IDE_4.Q(11,:),'LineWidth',1.5)
plot(IDE_4.Q(12,:),'LineWidth',1.5)
plot(IDE_4.Q(13,:),'LineWidth',1.5)
plot(IDE_4.Q(14,:),'LineWidth',1.5)
plot(IDE_4.Q(15,:),'LineWidth',1.5)
plot(IDE_4.Q(16,:),'LineWidth',1.5,'Color',[0 0 1])
xlabel('k-th trajectory point')
ylabel('joint value (m,rad)')
title('H) Trajectory 4')

Save_Figure(f,'Figure11')

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
f = figure;

subplot(4,2,1)
hold on
grid on
plot(KABC_1.Q(1,:),'LineWidth',1.5)
plot(KABC_1.Q(2,:),'LineWidth',1.5)
plot(KABC_1.Q(3,:),'LineWidth',1.5)
plot(KABC_1.Q(4,:),'LineWidth',1.5)
plot(KABC_1.Q(5,:),'LineWidth',1.5)
plot(KABC_1.Q(6,:),'LineWidth',1.5)
plot(KABC_1.Q(7,:),'LineWidth',1.5)
plot(KABC_1.Q(8,:),'LineWidth',1.5,'Color',[0 0 1])
xlabel('k-th trajectory point')
ylabel('joint value (m,rad)')
title('A) Trajectory 1')

subplot(4,2,2)
hold on
grid on
plot(KABC_2.Q(1,:),'LineWidth',1.5)
plot(KABC_2.Q(2,:),'LineWidth',1.5)
plot(KABC_2.Q(3,:),'LineWidth',1.5)
plot(KABC_2.Q(4,:),'LineWidth',1.5)
plot(KABC_2.Q(5,:),'LineWidth',1.5)
plot(KABC_2.Q(6,:),'LineWidth',1.5)
plot(KABC_2.Q(7,:),'LineWidth',1.5)
plot(KABC_2.Q(8,:),'LineWidth',1.5,'Color',[0 0 1])
xlabel('k-th trajectory point')
ylabel('joint value (m,rad)')
title('B) Trajectory 2')

subplot(4,2,3)
hold on
grid on
plot(KABC_3.Q(1,:),'LineWidth',1.5)
plot(KABC_3.Q(2,:),'LineWidth',1.5)
plot(KABC_3.Q(3,:),'LineWidth',1.5)
plot(KABC_3.Q(4,:),'LineWidth',1.5)
plot(KABC_3.Q(5,:),'LineWidth',1.5)
plot(KABC_3.Q(6,:),'LineWidth',1.5)
plot(KABC_3.Q(7,:),'LineWidth',1.5)
plot(KABC_3.Q(8,:),'LineWidth',1.5,'Color',[0 0 1])
xlabel('k-th trajectory point')
ylabel('joint value (m,rad)')
title('C) Trajectory 3')

subplot(4,2,4)
hold on
grid on
plot(KABC_4.Q(1,:),'LineWidth',1.5)
plot(KABC_4.Q(2,:),'LineWidth',1.5)
plot(KABC_4.Q(3,:),'LineWidth',1.5)
plot(KABC_4.Q(4,:),'LineWidth',1.5)
plot(KABC_4.Q(5,:),'LineWidth',1.5)
plot(KABC_4.Q(6,:),'LineWidth',1.5)
plot(KABC_4.Q(7,:),'LineWidth',1.5)
plot(KABC_4.Q(8,:),'LineWidth',1.5,'Color',[0 0 1])
xlabel('k-th trajectory point')
ylabel('joint value (m,rad)')
title('D) Trajectory 4')

subplot(4,2,5)
hold on
grid on
plot(KABC_1.Q(9,:),'LineWidth',1.5)
plot(KABC_1.Q(10,:),'LineWidth',1.5)
plot(KABC_1.Q(11,:),'LineWidth',1.5)
plot(KABC_1.Q(12,:),'LineWidth',1.5)
plot(KABC_1.Q(13,:),'LineWidth',1.5)
plot(KABC_1.Q(14,:),'LineWidth',1.5)
plot(KABC_1.Q(15,:),'LineWidth',1.5)
plot(KABC_1.Q(16,:),'LineWidth',1.5,'Color',[0 0 1])
xlabel('k-th trajectory point')
ylabel('joint value (m,rad)')
title('E) Trajectory 1')

subplot(4,2,6)
hold on
grid on
plot(KABC_2.Q(9,:),'LineWidth',1.5)
plot(KABC_2.Q(10,:),'LineWidth',1.5)
plot(KABC_2.Q(11,:),'LineWidth',1.5)
plot(KABC_2.Q(12,:),'LineWidth',1.5)
plot(KABC_2.Q(13,:),'LineWidth',1.5)
plot(KABC_2.Q(14,:),'LineWidth',1.5)
plot(KABC_2.Q(15,:),'LineWidth',1.5)
plot(KABC_2.Q(16,:),'LineWidth',1.5,'Color',[0 0 1])
xlabel('k-th trajectory point')
ylabel('joint value (m,rad)')
title('F) Trajectory 2')

subplot(4,2,7)
hold on
grid on
plot(KABC_3.Q(9,:),'LineWidth',1.5)
plot(KABC_3.Q(10,:),'LineWidth',1.5)
plot(KABC_3.Q(11,:),'LineWidth',1.5)
plot(KABC_3.Q(12,:),'LineWidth',1.5)
plot(KABC_3.Q(13,:),'LineWidth',1.5)
plot(KABC_3.Q(14,:),'LineWidth',1.5)
plot(KABC_3.Q(15,:),'LineWidth',1.5)
plot(KABC_3.Q(16,:),'LineWidth',1.5,'Color',[0 0 1])
xlabel('k-th trajectory point')
ylabel('joint value (m,rad)')
title('G) Trajectory 3')

subplot(4,2,8)
hold on
grid on
plot(KABC_4.Q(9,:),'LineWidth',1.5)
plot(KABC_4.Q(10,:),'LineWidth',1.5)
plot(KABC_4.Q(11,:),'LineWidth',1.5)
plot(KABC_4.Q(12,:),'LineWidth',1.5)
plot(KABC_4.Q(13,:),'LineWidth',1.5)
plot(KABC_4.Q(14,:),'LineWidth',1.5)
plot(KABC_4.Q(15,:),'LineWidth',1.5)
plot(KABC_4.Q(16,:),'LineWidth',1.5,'Color',[0 0 1])
xlabel('k-th trajectory point')
ylabel('joint value (m,rad)')
title('H) Trajectory 4')

Save_Figure(f,'Figure10')

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

function Save_Figure (f,name)
    ax = gcf;
    set (gcf, 'Position', [ax.Position(1)*0.1 ax.Position(2)*0.1 ax.Position(3) ax.Position(4)*2])

    ha = get(gcf,'children');
    
    set(ha(2),'position',[0.00+0.10 0.00+0.06 0.35 0.16]);
    set(ha(1),'position',[0.50+0.10 0.00+0.06 0.35 0.16]);
    
    set(ha(4),'position',[0.00+0.10 0.25+0.06 0.35 0.16]);
    set(ha(3),'position',[0.50+0.10 0.25+0.06 0.35 0.16]);

    set(ha(6),'position',[0.00+0.10 0.50+0.06 0.35 0.16]);
    set(ha(5),'position',[0.50+0.10 0.50+0.06 0.35 0.16]);

    set(ha(8),'position',[0.00+0.10 0.75+0.06 0.35 0.16]);
    set(ha(7),'position',[0.50+0.10 0.75+0.06 0.35 0.16]);

    legend(ha(7),{'$x_b$','$y_b$','$\theta_b$','$\theta_1$','$\theta_2$','$\theta_3$','$\theta_4$','$\theta_5$'},'Location','best','Interpreter','latex')
    
    print(f,['JointMotionResults/' name '.png'],'-dpng','-r300');
%     print(f,['JointMotionResults/' name '.eps'],'-depsc','-r300');
end
