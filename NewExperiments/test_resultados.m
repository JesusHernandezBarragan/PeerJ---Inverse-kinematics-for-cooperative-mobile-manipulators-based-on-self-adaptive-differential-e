clear all
close all
clc


S = 2;
T = 0.5; % default 0.5

f = figure;


load(['IDE_1_' num2str(S)])
[~,~,M] = size(Q);

for j=1:M
    subplot(2,2,1)
    hold on
    grid on
    plot(Q(1,:,j),'LineWidth',T,'Color',[0 0.4470 0.7410])
    plot(Q(2,:,j),'LineWidth',T,'Color',[0.8500 0.3250 0.0980])
    plot(Q(3,:,j),'LineWidth',T,'Color',[0.9290 0.6940 0.1250])
    plot(Q(4,:,j),'LineWidth',T,'Color',[0.4940 0.1840 0.5560])
    plot(Q(5,:,j),'LineWidth',T,'Color',[0.4660 0.6740 0.1880])
    plot(Q(6,:,j),'LineWidth',T,'Color',[0.3010 0.7450 0.9330])
    plot(Q(7,:,j),'LineWidth',T,'Color',[0.6350 0.0780 0.1840])
    plot(Q(8,:,j),'LineWidth',T,'Color',[0 0 1])
end

xlabel('k-th trajectory point')
ylabel('joint value (m,rad)')
title('Trajectory 1')


load(['IDE_2_' num2str(S)])
[~,~,M] = size(Q);

for j=1:M
    subplot(2,2,2)
    hold on
    grid on
    plot(Q(1,:,j),'LineWidth',T,'Color',[0 0.4470 0.7410])
    plot(Q(2,:,j),'LineWidth',T,'Color',[0.8500 0.3250 0.0980])
    plot(Q(3,:,j),'LineWidth',T,'Color',[0.9290 0.6940 0.1250])
    plot(Q(4,:,j),'LineWidth',T,'Color',[0.4940 0.1840 0.5560])
    plot(Q(5,:,j),'LineWidth',T,'Color',[0.4660 0.6740 0.1880])
    plot(Q(6,:,j),'LineWidth',T,'Color',[0.3010 0.7450 0.9330])
    plot(Q(7,:,j),'LineWidth',T,'Color',[0.6350 0.0780 0.1840])
    plot(Q(8,:,j),'LineWidth',T,'Color',[0 0 1])
    [~, hobj, ~, ~] = legend({'$x_b$','$y_b$','$\theta_b$','$\theta_1$','$\theta_2$','$\theta_3$','$\theta_4$','$\theta_5$'},'Interpreter','latex','Location','best');
end

hl = findobj(hobj,'type','line');
set(hl,'LineWidth',1.5);

xlabel('k-th trajectory point')
ylabel('joint value (m,rad)')
title('Trajectory 2')


load(['IDE_3_' num2str(S)])
[~,~,M] = size(Q);

for j=1:M
    subplot(2,2,3)
    hold on
    grid on
    plot(Q(1,:,j),'LineWidth',T,'Color',[0 0.4470 0.7410])
    plot(Q(2,:,j),'LineWidth',T,'Color',[0.8500 0.3250 0.0980])
    plot(Q(3,:,j),'LineWidth',T,'Color',[0.9290 0.6940 0.1250])
    plot(Q(4,:,j),'LineWidth',T,'Color',[0.4940 0.1840 0.5560])
    plot(Q(5,:,j),'LineWidth',T,'Color',[0.4660 0.6740 0.1880])
    plot(Q(6,:,j),'LineWidth',T,'Color',[0.3010 0.7450 0.9330])
    plot(Q(7,:,j),'LineWidth',T,'Color',[0.6350 0.0780 0.1840])
    plot(Q(8,:,j),'LineWidth',T,'Color',[0 0 1])
end

xlabel('k-th trajectory point')
ylabel('joint value (m,rad)')
title('Trajectory 3')



load(['IDE_4_' num2str(S)])
[~,~,M] = size(Q);

for j=1:M
    subplot(2,2,4)
    hold on
    grid on
    plot(Q(1,:,j),'LineWidth',T,'Color',[0 0.4470 0.7410])
    plot(Q(2,:,j),'LineWidth',T,'Color',[0.8500 0.3250 0.0980])
    plot(Q(3,:,j),'LineWidth',T,'Color',[0.9290 0.6940 0.1250])
    plot(Q(4,:,j),'LineWidth',T,'Color',[0.4940 0.1840 0.5560])
    plot(Q(5,:,j),'LineWidth',T,'Color',[0.4660 0.6740 0.1880])
    plot(Q(6,:,j),'LineWidth',T,'Color',[0.3010 0.7450 0.9330])
    plot(Q(7,:,j),'LineWidth',T,'Color',[0.6350 0.0780 0.1840])
    plot(Q(8,:,j),'LineWidth',T,'Color',[0 0 1])
end

xlabel('k-th trajectory point')
ylabel('joint value (m,rad)')
title('Trajectory 4')

Save_Figure(f,'IDE_MJointMotion_MM1')




f = figure;

load(['IDE_1_' num2str(S)])
[~,~,M] = size(Q);

for j=1:M
    subplot(2,2,1)
    hold on
    grid on
    plot(Q(9,:,j),'LineWidth',T,'Color',[0 0.4470 0.7410])
    plot(Q(10,:,j),'LineWidth',T,'Color',[0.8500 0.3250 0.0980])
    plot(Q(11,:,j),'LineWidth',T,'Color',[0.9290 0.6940 0.1250])
    plot(Q(12,:,j),'LineWidth',T,'Color',[0.4940 0.1840 0.5560])
    plot(Q(13,:,j),'LineWidth',T,'Color',[0.4660 0.6740 0.1880])
    plot(Q(14,:,j),'LineWidth',T,'Color',[0.3010 0.7450 0.9330])
    plot(Q(15,:,j),'LineWidth',T,'Color',[0.6350 0.0780 0.1840])
    plot(Q(16,:,j),'LineWidth',T,'Color',[0 0 1])
end

xlabel('k-th trajectory point')
ylabel('joint value (m,rad)')
title('Trajectory 1')


load(['IDE_2_' num2str(S)])
[~,~,M] = size(Q);

for j=1:M
    subplot(2,2,2)
    hold on
    grid on
    plot(Q(9,:,j),'LineWidth',T,'Color',[0 0.4470 0.7410])
    plot(Q(10,:,j),'LineWidth',T,'Color',[0.8500 0.3250 0.0980])
    plot(Q(11,:,j),'LineWidth',T,'Color',[0.9290 0.6940 0.1250])
    plot(Q(12,:,j),'LineWidth',T,'Color',[0.4940 0.1840 0.5560])
    plot(Q(13,:,j),'LineWidth',T,'Color',[0.4660 0.6740 0.1880])
    plot(Q(14,:,j),'LineWidth',T,'Color',[0.3010 0.7450 0.9330])
    plot(Q(15,:,j),'LineWidth',T,'Color',[0.6350 0.0780 0.1840])
    plot(Q(16,:,j),'LineWidth',T,'Color',[0 0 1])
    [~, hobj, ~, ~] = legend({'$x_b$','$y_b$','$\theta_b$','$\theta_1$','$\theta_2$','$\theta_3$','$\theta_4$','$\theta_5$'},'Interpreter','latex','Location','best');
end

hl = findobj(hobj,'type','line');
set(hl,'LineWidth',1.5);

xlabel('k-th trajectory point')
ylabel('joint value (m,rad)')
title('Trajectory 2')


load(['IDE_3_' num2str(S)])
[~,~,M] = size(Q);

for j=1:M
    subplot(2,2,3)
    hold on
    grid on
    plot(Q(9,:,j),'LineWidth',T,'Color',[0 0.4470 0.7410])
    plot(Q(10,:,j),'LineWidth',T,'Color',[0.8500 0.3250 0.0980])
    plot(Q(11,:,j),'LineWidth',T,'Color',[0.9290 0.6940 0.1250])
    plot(Q(12,:,j),'LineWidth',T,'Color',[0.4940 0.1840 0.5560])
    plot(Q(13,:,j),'LineWidth',T,'Color',[0.4660 0.6740 0.1880])
    plot(Q(14,:,j),'LineWidth',T,'Color',[0.3010 0.7450 0.9330])
    plot(Q(15,:,j),'LineWidth',T,'Color',[0.6350 0.0780 0.1840])
    plot(Q(16,:,j),'LineWidth',T,'Color',[0 0 1])
end

xlabel('k-th trajectory point')
ylabel('joint value (m,rad)')
title('Trajectory 3')



load(['IDE_4_' num2str(S)])
[~,~,M] = size(Q);

for j=1:M
    subplot(2,2,4)
    hold on
    grid on
    plot(Q(9,:,j),'LineWidth',T,'Color',[0 0.4470 0.7410])
    plot(Q(10,:,j),'LineWidth',T,'Color',[0.8500 0.3250 0.0980])
    plot(Q(11,:,j),'LineWidth',T,'Color',[0.9290 0.6940 0.1250])
    plot(Q(12,:,j),'LineWidth',T,'Color',[0.4940 0.1840 0.5560])
    plot(Q(13,:,j),'LineWidth',T,'Color',[0.4660 0.6740 0.1880])
    plot(Q(14,:,j),'LineWidth',T,'Color',[0.3010 0.7450 0.9330])
    plot(Q(15,:,j),'LineWidth',T,'Color',[0.6350 0.0780 0.1840])
    plot(Q(16,:,j),'LineWidth',T,'Color',[0 0 1])
end

xlabel('k-th trajectory point')
ylabel('joint value (m,rad)')
title('Trajectory 4')

Save_Figure(f,'IDE_MJointMotion_MM2')
