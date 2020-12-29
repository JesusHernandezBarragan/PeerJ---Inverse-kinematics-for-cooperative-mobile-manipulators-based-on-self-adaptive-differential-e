clear all
close all
clc

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

load IDE_4_2

[D,N,M] = size(Q);


M5 = zeros(N,M);


for i=1:N
    for j=1:M
        M5(i,j) = Q(5,i,j);
    end
end


figure
for j=1:M 

    hold on
    grid on
    plot(Q(1,:,j),'y-.')
    plot(Q(2,:,j),'m-.')
    plot(Q(3,:,j),'c-.')
    plot(Q(4,:,j),'r-.')
    plot(Q(5,:,j),'g-.')
    plot(Q(6,:,j),'b-.')
    plot(Q(7,:,j),'w-.')
    plot(Q(8,:,j),'k-.')
    legend('qa1','qa2','qa3','qa4','qa5','qa6','qa7','qa8')
end

figure
for j=1:M
    hold on
    grid on
    plot(Q(9,:,j),'y-.')
    plot(Q(10,:,j),'m-.')
    plot(Q(11,:,j),'c-.')
    plot(Q(12,:,j),'r-.')
    plot(Q(13,:,j),'g-.')
    plot(Q(14,:,j),'b-.')
    plot(Q(15,:,j),'w-.')
    plot(Q(16,:,j),'k-.')
    legend('qb1','qb2','qb3','qb4','qb5','qb6','qb7','qb8')
end



return

e = zeros(2,N);
td1 = zeros(3,N);
td2 = zeros(3,N);
t1 = zeros(3,N);
t2 = zeros(3,N);

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
end

figure
hold on
grid on
plot3(td1(1,:),td1(2,:),td1(3,:))
plot3(t1(1,:),t1(2,:),t1(3,:))
xlabel('x')
ylabel('y')
zlabel('z')
xlim([0.44 0.55])
view([-50 40])
legend('td1','t1')

figure
hold on
grid on
plot3(td2(1,:),td2(2,:),td2(3,:))
plot3(t2(1,:),t2(2,:),t2(3,:))
xlabel('x')
ylabel('y')
zlabel('z')
xlim([0.44 0.55])
view([-50 40])
legend('td2','t2')

figure
hold on
grid on
plot(e(1,:))
plot(e(2,:))
legend('e1','e2')



figure
boxplot(time)

[~,igb] = min(min(fitness));

figure
hold on
grid on
plot(fitness(:,igb))


