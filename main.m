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


% [mean([CFPSO_1.eq CFPSO_2.eq CFPSO_3.eq CFPSO_4.eq]);
% mean([FPA_1.eq FPA_2.eq FPA_3.eq FPA_4.eq]);
% mean([IDE_1.eq IDE_2.eq IDE_3.eq IDE_4.eq]);
% mean([KABC_1.eq KABC_2.eq KABC_3.eq KABC_4.eq])]

% return

% [mean([CFPSO_1.time CFPSO_2.time CFPSO_3.time CFPSO_4.time]);
% mean([FPA_1.time FPA_2.time FPA_3.time FPA_4.time]);
% mean([IDE_1.time IDE_2.time IDE_3.time IDE_4.time]);
% mean([KABC_1.time KABC_2.time KABC_3.time KABC_4.time])]

% return


% format shortE
% format compact

% [mean std min max] [t1_e1 t1_e2 t2_e1 t2_e2 t3_e1 t3_e2 t4_e1 t4_e2]
T = [[mean(CFPSO_1.e') mean(CFPSO_2.e') mean(CFPSO_3.e') mean(CFPSO_4.e')]; ...
[mean(FPA_1.e') mean(FPA_2.e') mean(FPA_3.e') mean(FPA_4.e')]; ...
[mean(IDE_1.e') mean(IDE_2.e') mean(IDE_3.e') mean(IDE_4.e')]; ...
[mean(KABC_1.e') mean(KABC_2.e') mean(KABC_3.e') mean(KABC_4.e')];];
T'

% [mean std min max] [t1_e1 t1_e2 t2_e1 t2_e2 t3_e1 t3_e2 t4_e1 t4_e2]
% T1 = [[mean(CFPSO_1.e')' std(CFPSO_1.e')' min(CFPSO_1.e')' max(CFPSO_1.e')']; ...
% [mean(FPA_1.e')' std(FPA_1.e')' min(FPA_1.e')' max(FPA_1.e')']; ...
% [mean(IDE_1.e')' std(IDE_1.e')' min(IDE_1.e')' max(IDE_1.e')']; ...
% [mean(KABC_1.e')' std(KABC_1.e')' min(KABC_1.e')' max(KABC_1.e')']];
% 
% T2 = [[mean(CFPSO_2.e')' std(CFPSO_2.e')' min(CFPSO_2.e')' max(CFPSO_2.e')']; ...
% [mean(FPA_2.e')' std(FPA_2.e')' min(FPA_2.e')' max(FPA_2.e')']; ...
% [mean(IDE_2.e')' std(IDE_2.e')' min(IDE_2.e')' max(IDE_2.e')']; ...
% [mean(KABC_2.e')' std(KABC_2.e')' min(KABC_2.e')' max(KABC_2.e')']];
% 
% T3 = [[mean(CFPSO_3.e')' std(CFPSO_3.e')' min(CFPSO_3.e')' max(CFPSO_3.e')']; ...
% [mean(FPA_3.e')' std(FPA_3.e')' min(FPA_3.e')' max(FPA_3.e')']; ...
% [mean(IDE_3.e')' std(IDE_3.e')' min(IDE_3.e')' max(IDE_3.e')']; ...
% [mean(KABC_3.e')' std(KABC_3.e')' min(KABC_3.e')' max(KABC_3.e')']];
% 
% T4 = [[mean(CFPSO_4.e')' std(CFPSO_4.e')' min(CFPSO_4.e')' max(CFPSO_4.e')']; ...
% [mean(FPA_4.e')' std(FPA_4.e')' min(FPA_4.e')' max(FPA_4.e')']; ...
% [mean(IDE_4.e')' std(IDE_4.e')' min(IDE_4.e')' max(IDE_4.e')']; ...
% [mean(KABC_4.e')' std(KABC_4.e')' min(KABC_4.e')' max(KABC_4.e')']];

format

return

%% Positon error
f = figure;

subplot(2,2,1)
boxplot([sum(CFPSO_1.e)' sum(FPA_1.e)' sum(IDE_1.e)' sum(KABC_1.e)'],'notch','off','labels',{s1,s2,s3,s4})
ylabel('error value (m)');
title('Trajectory 1');

subplot(2,2,2)
boxplot([sum(CFPSO_2.e)' sum(FPA_2.e)' sum(IDE_2.e)' sum(KABC_2.e)'],'notch','off','labels',{s1,s2,s3,s4})
ylabel('error value (m)');
title('Trajectory 2');

subplot(2,2,3)
boxplot([sum(CFPSO_3.e)' sum(FPA_3.e)' sum(IDE_3.e)' sum(KABC_3.e)'],'notch','off','labels',{s1,s2,s3,s4})
ylabel('error value (m)');
title('Trajectory 3');

subplot(2,2,4)
boxplot([sum(CFPSO_4.e)' sum(FPA_4.e)' sum(IDE_4.e)' sum(KABC_4.e)'],'notch','off','labels',{s1,s2,s3,s4})
ylabel('error value (m)');
title('Trajectory 4');

Save_Figure(f,'PositionError')

%% Motion error
f = figure;

subplot(2,2,1)
boxplot([CFPSO_1.eq' FPA_1.eq' IDE_1.eq' KABC_1.eq'],'notch','off','labels',{s1,s2,s3,s4})
ylabel('error value');
title('Trajectory 1');

subplot(2,2,2)
boxplot([CFPSO_2.eq' FPA_2.eq' IDE_2.eq' KABC_2.eq'],'notch','off','labels',{s1,s2,s3,s4})
ylabel('error value');
title('Trajectory 2');

subplot(2,2,3)
boxplot([CFPSO_3.eq' FPA_3.eq' IDE_3.eq' KABC_3.eq'],'notch','off','labels',{s1,s2,s3,s4})
ylabel('error value');
title('Trajectory 3');

subplot(2,2,4)
boxplot([CFPSO_4.eq' FPA_4.eq' IDE_4.eq' KABC_4.eq'],'notch','off','labels',{s1,s2,s3,s4})
ylabel('error value');
title('Trajectory 4');

Save_Figure(f,'MotionError')

%% Execution Time
f = figure;

subplot(2,2,1)
boxplot([CFPSO_1.time' FPA_1.time' IDE_1.time' KABC_1.time'],'notch','off','labels',{s1,s2,s3,s4})
ylabel('time (s)');
title('Trajectory 1');

subplot(2,2,2)
boxplot([CFPSO_2.time' FPA_2.time' IDE_2.time' KABC_2.time'],'notch','off','labels',{s1,s2,s3,s4})
ylabel('time (s)');
title('Trajectory 2');

subplot(2,2,3)
boxplot([CFPSO_3.time' FPA_3.time' IDE_3.time' KABC_3.time'],'notch','off','labels',{s1,s2,s3,s4})
ylabel('time (s)');
title('Trajectory 3');

subplot(2,2,4)
boxplot([CFPSO_4.time' FPA_4.time' IDE_4.time' KABC_4.time'],'notch','off','labels',{s1,s2,s3,s4})
ylabel('time (s)');
title('Trajectory 4');

Save_Figure(f,'ExecutionTime')

%% Position error per iteration
% r1 = randi([1 CFPSO_1.N]);
% r2 = randi([1 FPA_1.N]);
% r3 = randi([1 IDE_1.N]);
% r4 = randi([1 KABC_1.N]);
% yup = 0.4;
% 
% f = figure;
% 
% subplot(2,2,1)
% hold on
% grid on
% plot(CFPSO_1.fit(:,r1),'-','LineWidth',1.5)
% plot(FPA_1.fit(:,r2),'--','LineWidth',1.5)
% plot(IDE_1.fit(:,r3),'-.','LineWidth',1.5)
% plot(KABC_1.fit(:,r4),':','LineWidth',1.5)
% title('Trajectory 1');
% xlabel('iteration')
% ylabel('position error (m)');
% % legend(s1,s2,s3,s4)
% ylim([0 yup])
% 
% subplot(2,2,2)
% hold on
% grid on
% plot(CFPSO_2.fit(:,r1),'-','LineWidth',1.5)
% plot(FPA_2.fit(:,r2),'--','LineWidth',1.5)
% plot(IDE_2.fit(:,r3),'-.','LineWidth',1.5)
% plot(KABC_2.fit(:,r4),':','LineWidth',1.5)
% title('Trajectory 2');
% xlabel('iteration')
% ylabel('position error (m)');
% legend(s1,s2,s3,s4)
% ylim([0 yup])
% 
% subplot(2,2,3)
% hold on
% grid on
% plot(CFPSO_3.fit(:,r1),'-','LineWidth',1.5)
% plot(FPA_3.fit(:,r2),'--','LineWidth',1.5)
% plot(IDE_3.fit(:,r3),'-.','LineWidth',1.5)
% plot(KABC_3.fit(:,r4),':','LineWidth',1.5)
% title('Trajectory 3');
% xlabel('iteration')
% ylabel('position error (m)');
% % legend(s1,s2,s3,s4)
% ylim([0 yup])
% 
% subplot(2,2,4)
% hold on
% grid on
% plot(CFPSO_4.fit(:,r1),'-','LineWidth',1.5)
% plot(FPA_4.fit(:,r2),'--','LineWidth',1.5)
% plot(IDE_4.fit(:,r3),'-.','LineWidth',1.5)
% plot(KABC_4.fit(:,r4),':','LineWidth',1.5)
% title('Trajectory 4');
% xlabel('iteration')
% ylabel('position error (m)');
% % legend(s1,s2,s3,s4)
% ylim([0 yup])
% 
% Save_Figure(f,'IterationError')

%% Joint Motion 
f = figure;

subplot(2,2,1)
hold on
grid on
plot(IDE_1.Q(1,:),'LineWidth',1.5)
plot(IDE_1.Q(2,:),'LineWidth',1.5)
plot(IDE_1.Q(3,:),'LineWidth',1.5)
plot(IDE_1.Q(4,:),'LineWidth',1.5)
plot(IDE_1.Q(5,:),'LineWidth',1.5)
plot(IDE_1.Q(6,:),'LineWidth',1.5)
plot(IDE_1.Q(7,:),'LineWidth',1.5)
plot(IDE_1.Q(8,:),'LineWidth',1.5)
xlabel('k-th trajectory point')
ylabel('joint value (m,rad)')
title('Trajectory 1')

subplot(2,2,2)
hold on
grid on
plot(IDE_2.Q(1,:),'LineWidth',1.5)
plot(IDE_2.Q(2,:),'LineWidth',1.5)
plot(IDE_2.Q(3,:),'LineWidth',1.5)
plot(IDE_2.Q(4,:),'LineWidth',1.5)
plot(IDE_2.Q(5,:),'LineWidth',1.5)
plot(IDE_2.Q(6,:),'LineWidth',1.5)
plot(IDE_2.Q(7,:),'LineWidth',1.5)
plot(IDE_2.Q(8,:),'LineWidth',1.5)
legend({'$x_b$','$y_b$','$\theta_b$','$\theta_1$','$\theta_2$','$\theta_3$','$\theta_4$','$\theta_5$'},'Interpreter','latex','Location','best')
xlabel('k-th trajectory point')
ylabel('joint value (m,rad)')
title('Trajectory 2')

subplot(2,2,3)
hold on
grid on
plot(IDE_3.Q(1,:),'LineWidth',1.5)
plot(IDE_3.Q(2,:),'LineWidth',1.5)
plot(IDE_3.Q(3,:),'LineWidth',1.5)
plot(IDE_3.Q(4,:),'LineWidth',1.5)
plot(IDE_3.Q(5,:),'LineWidth',1.5)
plot(IDE_3.Q(6,:),'LineWidth',1.5)
plot(IDE_3.Q(7,:),'LineWidth',1.5)
plot(IDE_3.Q(8,:),'LineWidth',1.5)
xlabel('k-th trajectory point')
ylabel('joint value (m,rad)')
title('Trajectory 3')

subplot(2,2,4)
hold on
grid on
plot(IDE_4.Q(1,:),'LineWidth',1.5)
plot(IDE_4.Q(2,:),'LineWidth',1.5)
plot(IDE_4.Q(3,:),'LineWidth',1.5)
plot(IDE_4.Q(4,:),'LineWidth',1.5)
plot(IDE_4.Q(5,:),'LineWidth',1.5)
plot(IDE_4.Q(6,:),'LineWidth',1.5)
plot(IDE_4.Q(7,:),'LineWidth',1.5)
plot(IDE_4.Q(8,:),'LineWidth',1.5)
xlabel('k-th trajectory point')
ylabel('joint value (m,rad)')
title('Trajectory 4')

Save_Figure(f,'IDE_JointMotion_MM1')


f = figure;

subplot(2,2,1)
hold on
grid on
plot(IDE_1.Q(9,:),'LineWidth',1.5)
plot(IDE_1.Q(10,:),'LineWidth',1.5)
plot(IDE_1.Q(11,:),'LineWidth',1.5)
plot(IDE_1.Q(12,:),'LineWidth',1.5)
plot(IDE_1.Q(13,:),'LineWidth',1.5)
plot(IDE_1.Q(14,:),'LineWidth',1.5)
plot(IDE_1.Q(15,:),'LineWidth',1.5)
plot(IDE_1.Q(16,:),'LineWidth',1.5)
xlabel('k-th trajectory point')
ylabel('joint value (m,rad)')
title('Trajectory 1')

subplot(2,2,2)
hold on
grid on
plot(IDE_2.Q(9,:),'LineWidth',1.5)
plot(IDE_2.Q(10,:),'LineWidth',1.5)
plot(IDE_2.Q(11,:),'LineWidth',1.5)
plot(IDE_2.Q(12,:),'LineWidth',1.5)
plot(IDE_2.Q(13,:),'LineWidth',1.5)
plot(IDE_2.Q(14,:),'LineWidth',1.5)
plot(IDE_2.Q(15,:),'LineWidth',1.5)
plot(IDE_2.Q(16,:),'LineWidth',1.5)
legend({'$x_b$','$y_b$','$\theta_b$','$\theta_1$','$\theta_2$','$\theta_3$','$\theta_4$','$\theta_5$'},'Interpreter','latex','Location','best')
xlabel('k-th trajectory point')
ylabel('joint value (m,rad)')
title('Trajectory 2')

subplot(2,2,3)
hold on
grid on
plot(IDE_3.Q(9,:),'LineWidth',1.5)
plot(IDE_3.Q(10,:),'LineWidth',1.5)
plot(IDE_3.Q(11,:),'LineWidth',1.5)
plot(IDE_3.Q(12,:),'LineWidth',1.5)
plot(IDE_3.Q(13,:),'LineWidth',1.5)
plot(IDE_3.Q(14,:),'LineWidth',1.5)
plot(IDE_3.Q(15,:),'LineWidth',1.5)
plot(IDE_3.Q(16,:),'LineWidth',1.5)
xlabel('k-th trajectory point')
ylabel('joint value (m,rad)')
title('Trajectory 3')

subplot(2,2,4)
hold on
grid on
plot(IDE_4.Q(9,:),'LineWidth',1.5)
plot(IDE_4.Q(10,:),'LineWidth',1.5)
plot(IDE_4.Q(11,:),'LineWidth',1.5)
plot(IDE_4.Q(12,:),'LineWidth',1.5)
plot(IDE_4.Q(13,:),'LineWidth',1.5)
plot(IDE_4.Q(14,:),'LineWidth',1.5)
plot(IDE_4.Q(15,:),'LineWidth',1.5)
plot(IDE_4.Q(16,:),'LineWidth',1.5)
xlabel('k-th trajectory point')
ylabel('joint value (m,rad)')
title('Trajectory 4')

Save_Figure(f,'IDE_JointMotion_MM2')

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
f = figure;

subplot(2,2,1)
hold on
grid on
plot(KABC_1.Q(1,:),'LineWidth',1.5)
plot(KABC_1.Q(2,:),'LineWidth',1.5)
plot(KABC_1.Q(3,:),'LineWidth',1.5)
plot(KABC_1.Q(4,:),'LineWidth',1.5)
plot(KABC_1.Q(5,:),'LineWidth',1.5)
plot(KABC_1.Q(6,:),'LineWidth',1.5)
plot(KABC_1.Q(7,:),'LineWidth',1.5)
plot(KABC_1.Q(8,:),'LineWidth',1.5)
xlabel('k-th trajectory point')
ylabel('joint value (m,rad)')
title('Trajectory 1')

subplot(2,2,2)
hold on
grid on
plot(KABC_2.Q(1,:),'LineWidth',1.5)
plot(KABC_2.Q(2,:),'LineWidth',1.5)
plot(KABC_2.Q(3,:),'LineWidth',1.5)
plot(KABC_2.Q(4,:),'LineWidth',1.5)
plot(KABC_2.Q(5,:),'LineWidth',1.5)
plot(KABC_2.Q(6,:),'LineWidth',1.5)
plot(KABC_2.Q(7,:),'LineWidth',1.5)
plot(KABC_2.Q(8,:),'LineWidth',1.5)
legend({'$x_b$','$y_b$','$\theta_b$','$\theta_1$','$\theta_2$','$\theta_3$','$\theta_4$','$\theta_5$'},'Interpreter','latex','Location','best')
xlabel('k-th trajectory point')
ylabel('joint value (m,rad)')
title('Trajectory 2')

subplot(2,2,3)
hold on
grid on
plot(KABC_3.Q(1,:),'LineWidth',1.5)
plot(KABC_3.Q(2,:),'LineWidth',1.5)
plot(KABC_3.Q(3,:),'LineWidth',1.5)
plot(KABC_3.Q(4,:),'LineWidth',1.5)
plot(KABC_3.Q(5,:),'LineWidth',1.5)
plot(KABC_3.Q(6,:),'LineWidth',1.5)
plot(KABC_3.Q(7,:),'LineWidth',1.5)
plot(KABC_3.Q(8,:),'LineWidth',1.5)
xlabel('k-th trajectory point')
ylabel('joint value (m,rad)')
title('Trajectory 3')

subplot(2,2,4)
hold on
grid on
plot(KABC_4.Q(1,:),'LineWidth',1.5)
plot(KABC_4.Q(2,:),'LineWidth',1.5)
plot(KABC_4.Q(3,:),'LineWidth',1.5)
plot(KABC_4.Q(4,:),'LineWidth',1.5)
plot(KABC_4.Q(5,:),'LineWidth',1.5)
plot(KABC_4.Q(6,:),'LineWidth',1.5)
plot(KABC_4.Q(7,:),'LineWidth',1.5)
plot(KABC_4.Q(8,:),'LineWidth',1.5)
xlabel('k-th trajectory point')
ylabel('joint value (m,rad)')
title('Trajectory 4')

Save_Figure(f,'KABC_JointMotion_MM1')


f = figure;

subplot(2,2,1)
hold on
grid on
plot(KABC_1.Q(9,:),'LineWidth',1.5)
plot(KABC_1.Q(10,:),'LineWidth',1.5)
plot(KABC_1.Q(11,:),'LineWidth',1.5)
plot(KABC_1.Q(12,:),'LineWidth',1.5)
plot(KABC_1.Q(13,:),'LineWidth',1.5)
plot(KABC_1.Q(14,:),'LineWidth',1.5)
plot(KABC_1.Q(15,:),'LineWidth',1.5)
plot(KABC_1.Q(16,:),'LineWidth',1.5)
xlabel('k-th trajectory point')
ylabel('joint value (m,rad)')
title('Trajectory 1')

subplot(2,2,2)
hold on
grid on
plot(KABC_2.Q(9,:),'LineWidth',1.5)
plot(KABC_2.Q(10,:),'LineWidth',1.5)
plot(KABC_2.Q(11,:),'LineWidth',1.5)
plot(KABC_2.Q(12,:),'LineWidth',1.5)
plot(KABC_2.Q(13,:),'LineWidth',1.5)
plot(KABC_2.Q(14,:),'LineWidth',1.5)
plot(KABC_2.Q(15,:),'LineWidth',1.5)
plot(KABC_2.Q(16,:),'LineWidth',1.5)
legend({'$x_b$','$y_b$','$\theta_b$','$\theta_1$','$\theta_2$','$\theta_3$','$\theta_4$','$\theta_5$'},'Interpreter','latex','Location','best')
xlabel('k-th trajectory point')
ylabel('joint value (m,rad)')
title('Trajectory 2')

subplot(2,2,3)
hold on
grid on
plot(KABC_3.Q(9,:),'LineWidth',1.5)
plot(KABC_3.Q(10,:),'LineWidth',1.5)
plot(KABC_3.Q(11,:),'LineWidth',1.5)
plot(KABC_3.Q(12,:),'LineWidth',1.5)
plot(KABC_3.Q(13,:),'LineWidth',1.5)
plot(KABC_3.Q(14,:),'LineWidth',1.5)
plot(KABC_3.Q(15,:),'LineWidth',1.5)
plot(KABC_3.Q(16,:),'LineWidth',1.5)
xlabel('k-th trajectory point')
ylabel('joint value (m,rad)')
title('Trajectory 3')

subplot(2,2,4)
hold on
grid on
plot(KABC_4.Q(9,:),'LineWidth',1.5)
plot(KABC_4.Q(10,:),'LineWidth',1.5)
plot(KABC_4.Q(11,:),'LineWidth',1.5)
plot(KABC_4.Q(12,:),'LineWidth',1.5)
plot(KABC_4.Q(13,:),'LineWidth',1.5)
plot(KABC_4.Q(14,:),'LineWidth',1.5)
plot(KABC_4.Q(15,:),'LineWidth',1.5)
plot(KABC_4.Q(16,:),'LineWidth',1.5)
xlabel('k-th trajectory point')
ylabel('joint value (m,rad)')
title('Trajectory 4')

Save_Figure(f,'KABC_JointMotion_MM2')

%% Tracking results
f = figure;

subplot(2,2,1)
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
title('Trajectory 1')

subplot(2,2,2)
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
title('Trajectory 2')
legend({'$\mathbf{t}_{e_{1}}^{*}$','$\mathbf{t}_{e_{1}}^{1}$'},'Interpreter','latex','Location','best','Location','best')

subplot(2,2,3)
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
title('Trajectory 3')

subplot(2,2,4)
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
title('Trajectory 4')

Save_Figure(f,'IDE_Tracking_MM1')


f = figure;

subplot(2,2,1)
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
title('Trajectory 1')

subplot(2,2,2)
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
title('Trajectory 2')
legend({'$\mathbf{t}_{e_{2}}^{*}$','$\mathbf{t}_{e_{2}}^{1}$'},'Interpreter','latex','Location','best')

subplot(2,2,3)
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
title('Trajectory 3')

subplot(2,2,4)
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
title('Trajectory 4')

Save_Figure(f,'IDE_Tracking_MM2')




% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

f = figure;

subplot(2,2,1)
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
title('Trajectory 1')

subplot(2,2,2)
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
title('Trajectory 2')
legend({'$\mathbf{t}_{e_{1}}^{*}$','$\mathbf{t}_{e_{1}}^{1}$'},'Interpreter','latex','Location','best')

subplot(2,2,3)
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
title('Trajectory 3')

subplot(2,2,4)
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
title('Trajectory 4')

Save_Figure(f,'KABC_Tracking_MM1')


f = figure;

subplot(2,2,1)
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
title('Trajectory 1')

subplot(2,2,2)
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
title('Trajectory 2')
legend({'$\mathbf{t}_{e_{2}}^{*}$','$\mathbf{t}_{e_{2}}^{1}$'},'Interpreter','latex','Location','best')

subplot(2,2,3)
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
title('Trajectory 3')

subplot(2,2,4)
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
title('Trajectory 4')

Save_Figure(f,'KABC_Tracking_MM2')

%%






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

% z1 = Sign_Test(sum(IDE_1.e),sum(CFPSO_1.e));
% z2 = Sign_Test(sum(IDE_1.e),sum(FPA_1.e));
% z3 = Sign_Test(sum(IDE_1.e),sum(KABC_1.e));
% Z1 = [z1 z2 z3];
% 
% z1 = Sign_Test(sum(IDE_2.e),sum(CFPSO_2.e));
% z2 = Sign_Test(sum(IDE_2.e),sum(FPA_2.e));
% z3 = Sign_Test(sum(IDE_2.e),sum(KABC_2.e));
% Z2 = [z1 z2 z3];
% 
% z1 = Sign_Test(sum(IDE_3.e),sum(CFPSO_3.e));
% z2 = Sign_Test(sum(IDE_3.e),sum(FPA_3.e));
% z3 = Sign_Test(sum(IDE_3.e),sum(KABC_3.e));
% Z3 = [z1 z2 z3];
% 
% z1 = Sign_Test(sum(IDE_4.e),sum(CFPSO_4.e));
% z2 = Sign_Test(sum(IDE_4.e),sum(FPA_4.e));
% z3 = Sign_Test(sum(IDE_4.e),sum(KABC_4.e));
% Z4 = [z1 z2 z3];

function z = Sign_Test (X,Y)
    n = numel(X);
    [p,~,stats] = signtest(X,Y);
    z = [n-stats.sign stats.sign p]';
end