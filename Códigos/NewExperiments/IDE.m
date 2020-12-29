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

prd = [0.0 -1.0 0.0]';

M = 25; % 25
Traj = 1;

N = 200;
[x,y,z] = Trayectoria(Traj,N);
Q = zeros(16,N,M);
time = zeros(1,N,M);
fitness = zeros(1000,N,M);

for j=1:M
    q0 = [+0.0 +0.5 0 0 +pi/2 -pi/4 +pi/4 0.0 ...
          +0.0 -0.5 0 0 +pi/2 -pi/4 +pi/4 0.0]'; 

%     q0 = [+0.3 +0.5 +pi/4 0 +pi/2 0 +pi/2 +pi/4 ...
%           +0.3 -0.5 -pi/4 0 +pi/2 0 +pi/2 -pi/4]'; 
  
    for i=1:N
        tic

        td1 = [x(i) y(i) z(i)]';
        td2 = td1 + prd;

        [xgb,fitness(:,i,j)] = IK (q0,td1,td2,Txy,Tyb,Tb0,T01,T12,T23,T34,T45);
        q0 = xgb;

        Q(:,i,j) = xgb;

        time(i,j) = toc;
        disp(['iter = ' num2str(i) ', time = ' num2str(time(i,j)) ', m = ' num2str(j)])
    end
end

save(['IDE_' num2str(Traj)],'x','y','z','prd','Q','N','M','Traj','time','fitness')



%% Funciones
function [xgb,plot_fit] = IK (q0,td1,td2,Txy,Tyb,Tb0,T01,T12,T23,T34,T45)
    tol = 0.0001;

    G = 1000;
    N = 30;
    D = 16;

    F = 0.6;
    CR = 0.9;

    xl = [-1.5 -1.5 -pi -169*(pi/180) -65*(pi/180) -150*(pi/180) -102.2*(pi/180) -167.5*(pi/180) -1.5 -1.5 -pi -169*(pi/180) -65*(pi/180) -150*(pi/180) -102.2*(pi/180) -167.5*(pi/180)]';
    xu = [+1.5 +1.5 +pi +169*(pi/180) +90*(pi/180) +146*(pi/180) +102.5*(pi/180) +167.5*(pi/180) +1.5 +1.5 +pi +169*(pi/180) +90*(pi/180) +146*(pi/180) +102.5*(pi/180) +167.5*(pi/180)]';

    x = zeros(D,N);
    fitness = zeros(1,N);

    f = @(t1,t2,q) 1.5*(norm(td1-t1)+norm(td2-t2)) + 0.6*norm(q0-q) + 100*Penalty(q,xl,xu);
    plot_fit = zeros(1,G);
    
    for i=1:N
        x(:,i) = xl+(xu-xl).*rand(D,1);

        q1 = x(1:8,i);
        Tx = Txy(q1(1));
        Ty = Tx*Tyb(q1(2));
        Tb = Ty*Tb0(q1(3));
        T1 = Tb*T01(q1(4));
        T2 = T1*T12(q1(5));
        T3 = T2*T23(q1(6));
        T4 = T3*T34(q1(7));
        T5 = T4*T45(q1(8));
        Tq1 = T5;

        q2 = x(9:16,i);
        Tx = Txy(q2(1));
        Ty = Tx*Tyb(q2(2));
        Tb = Ty*Tb0(q2(3));
        T1 = Tb*T01(q2(4));
        T2 = T1*T12(q2(5));
        T3 = T2*T23(q2(6));
        T4 = T3*T34(q2(7));
        T5 = T4*T45(q2(8));
        Tq2 = T5;

        fitness(i) = f(Tq1(1:3,4),Tq2(1:3,4),x(:,i));            
    end

    for g=1:G
        for i=1:N
            r1=i;
            while r1==i
                r1 = randi([1,N]);  
            end

            r2 = r1;
            while r2==r1 || r2==i
                r2 = randi([1,N]);
            end

            r3 = r2;
            while r3==r2 || r3==r1 || r3==i
                r3 = randi([1,N]);
            end

            r4 = r3;
            while r4==r3 || r4==r2 || r4==r1 || r4==i
                r4 = randi([1,N]);
            end
            
            if rand<(1-(i/N)^2)
                v = x(:,r1) + F*(x(:,r2)-x(:,r3));
            else
                [~,igb] = min(fitness);
                v = x(:,igb) + F*(x(:,r1)-x(:,r2)+x(:,r3)-x(:,r4));
            end
            
            u = zeros(D,1);
            j0 = randi([1 D]);

            for j=1:D
                if rand<=CR || j==j0
                    u(j) = v(j);
                else
                    u(j) = x(j,i);
                end
            end   

            q1 = u(1:8);
            Tx = Txy(q1(1));
            Ty = Tx*Tyb(q1(2));
            Tb = Ty*Tb0(q1(3));
            T1 = Tb*T01(q1(4));
            T2 = T1*T12(q1(5));
            T3 = T2*T23(q1(6));
            T4 = T3*T34(q1(7));
            T5 = T4*T45(q1(8));
            Tq1 = T5;

            q2 = u(9:16);
            Tx = Txy(q2(1));
            Ty = Tx*Tyb(q2(2));
            Tb = Ty*Tb0(q2(3));
            T1 = Tb*T01(q2(4));
            T2 = T1*T12(q2(5));
            T3 = T2*T23(q2(6));
            T4 = T3*T34(q2(7));
            T5 = T4*T45(q2(8));
            Tq2 = T5;

            fitness_u = f(Tq1(1:3,4),Tq2(1:3,4),u);

            if fitness_u < fitness(i)
                x(:,i) = u;
                fitness(i) = fitness_u;
            end  
        end

        [~,igb] = min(fitness);

        q1 = x(1:8,igb);
        Tx = Txy(q1(1));
        Ty = Tx*Tyb(q1(2));
        Tb = Ty*Tb0(q1(3));
        T1 = Tb*T01(q1(4));
        T2 = T1*T12(q1(5));
        T3 = T2*T23(q1(6));
        T4 = T3*T34(q1(7));
        T5 = T4*T45(q1(8));
        Tq1 = T5;

        q2 = x(9:16,igb);
        Tx = Txy(q2(1));
        Ty = Tx*Tyb(q2(2));
        Tb = Ty*Tb0(q2(3));
        T1 = Tb*T01(q2(4));
        T2 = T1*T12(q2(5));
        T3 = T2*T23(q2(6));
        T4 = T3*T34(q2(7));
        T5 = T4*T45(q2(8));
        Tq2 = T5;

        e1 = norm(td1-Tq1(1:3,4));
        e2 = norm(td2-Tq2(1:3,4));

        if(e1<tol && e2<tol)
            break
        end
        
%         plot_fit(g) = fitness(igb);
        plot_fit(g) = e1 + e2;
    end

    [~,igb] = min(fitness);
    xgb = x(:,igb);
end

function Dibujar_MM (q,Tx,Ty,Tb,T1,T2,T3,T4,T5,s)
    plot3(Tb(1,4),Tb(2,4),Tb(3,4),'bo','LineWidth',2,'MarkerSize',10)
    plot3(T1(1,4),T1(2,4),T1(3,4),'bo','LineWidth',2,'MarkerSize',10)
    plot3(T2(1,4),T2(2,4),T2(3,4),'bo','LineWidth',2,'MarkerSize',10)
    plot3(T3(1,4),T3(2,4),T3(3,4),'bo','LineWidth',2,'MarkerSize',10)
    plot3(T4(1,4),T4(2,4),T4(3,4),'bo','LineWidth',2,'MarkerSize',10)

    line([Ty(1,4) Tb(1,4)],[Ty(2,4) Tb(2,4)],[Ty(3,4) Tb(3,4)],'color',[0 0 1],'LineWidth',3,'MarkerSize',15)
    line([Tb(1,4) T1(1,4)],[Tb(2,4) T1(2,4)],[Tb(3,4) T1(3,4)],'color',[0 0 1],'LineWidth',3,'MarkerSize',15)
    line([T1(1,4) T2(1,4)],[T1(2,4) T2(2,4)],[T1(3,4) T2(3,4)],'color',[0 0 1],'LineWidth',3,'MarkerSize',15)
    line([T2(1,4) T3(1,4)],[T2(2,4) T3(2,4)],[T2(3,4) T3(3,4)],'color',[0 0 1],'LineWidth',3,'MarkerSize',15)
    line([T3(1,4) T4(1,4)],[T3(2,4) T4(2,4)],[T3(3,4) T4(3,4)],'color',[0 0 1],'LineWidth',3,'MarkerSize',15)
    line([T4(1,4) T5(1,4)],[T4(2,4) T5(2,4)],[T4(3,4) T5(3,4)],'color',[0 0 1],'LineWidth',3,'MarkerSize',15)

    Dibujar_Omnidireccional_4(q(1:3),0.25,0.20)
    Dibujar_Sistema_Referencia_3D(T5,s)
end

function Dibujar_Omnidireccional_4 (p,L,l)
    x = p(1);
    y = p(2);
    theta = p(3);

    hold on
    grid on
    xlabel('x')
    ylabel('y')
    
	Lo = L*0.3;
    lo = L*0.2;
    Tob = [cos(theta) -sin(theta) x; sin(theta) cos(theta) y; 0 0 1];
    Tblf = [1 0 L; 0 1 l; 0 0 1];
    Tbrf = [1 0 L; 0 1 -l; 0 0 1];
    Tblb = [1 0 -L; 0 1 l; 0 0 1];
    Tbrb = [1 0 -L; 0 1 -l; 0 0 1];
    
    Tolf = Tob*Tblf;
    Torf = Tob*Tbrf;
    Tolb = Tob*Tblb;
    Torb = Tob*Tbrb;
    
    % Base
    phi = linspace(0,2*pi,50);
    pc = Tob*[L*0.7 0 1]';
    
    cx = pc(1)+L*0.15*cos(phi); 
    cy = pc(2)+L*0.15*sin(phi);
    plot(cx,cy,'LineWidth',2,'MarkerSize',10,'color',[1 0 0])
    
    
    p1 = Tob*[+L+Lo -l+lo 1]';
    p2 = Tob*[-L-Lo -l+lo 1]';
    p3 = Tob*[+L+Lo +l-lo 1]';
    p4 = Tob*[-L-Lo +l-lo 1]';
    line([p1(1) p2(1)],[p1(2) p2(2)],'LineWidth',2,'MarkerSize',10,'color',[0 0 1])
    line([p1(1) p3(1)],[p1(2) p3(2)],'LineWidth',2,'MarkerSize',10,'color',[0 0 1])
    line([p2(1) p4(1)],[p2(2) p4(2)],'LineWidth',2,'MarkerSize',10,'color',[0 0 1])
    line([p3(1) p4(1)],[p3(2) p4(2)],'LineWidth',2,'MarkerSize',10,'color',[0 0 1])
    
    % Ruedas
    p1 = Tolf*[+Lo -lo 1]';
    p2 = Tolf*[-Lo -lo 1]';
    p3 = Tolf*[+Lo +lo 1]';
    p4 = Tolf*[-Lo +lo 1]';
    line([p1(1) p2(1)],[p1(2) p2(2)],'LineWidth',2,'MarkerSize',10,'color',[0 0 1])
    line([p1(1) p3(1)],[p1(2) p3(2)],'LineWidth',2,'MarkerSize',10,'color',[0 0 1])
    line([p2(1) p4(1)],[p2(2) p4(2)],'LineWidth',2,'MarkerSize',10,'color',[0 0 1])
    line([p3(1) p4(1)],[p3(2) p4(2)],'LineWidth',2,'MarkerSize',10,'color',[0 0 1])

    p1 = Torf*[+Lo -lo 1]';
    p2 = Torf*[-Lo -lo 1]';
    p3 = Torf*[+Lo +lo 1]';
    p4 = Torf*[-Lo +lo 1]';
    line([p1(1) p2(1)],[p1(2) p2(2)],'LineWidth',2,'MarkerSize',10,'color',[0 0 1])
    line([p1(1) p3(1)],[p1(2) p3(2)],'LineWidth',2,'MarkerSize',10,'color',[0 0 1])
    line([p2(1) p4(1)],[p2(2) p4(2)],'LineWidth',2,'MarkerSize',10,'color',[0 0 1])
    line([p3(1) p4(1)],[p3(2) p4(2)],'LineWidth',2,'MarkerSize',10,'color',[0 0 1])

    p1 = Torb*[+Lo -lo 1]';
    p2 = Torb*[-Lo -lo 1]';
    p3 = Torb*[+Lo +lo 1]';
    p4 = Torb*[-Lo +lo 1]';
    line([p1(1) p2(1)],[p1(2) p2(2)],'LineWidth',2,'MarkerSize',10,'color',[0 0 1])
    line([p1(1) p3(1)],[p1(2) p3(2)],'LineWidth',2,'MarkerSize',10,'color',[0 0 1])
    line([p2(1) p4(1)],[p2(2) p4(2)],'LineWidth',2,'MarkerSize',10,'color',[0 0 1])
    line([p3(1) p4(1)],[p3(2) p4(2)],'LineWidth',2,'MarkerSize',10,'color',[0 0 1])
    
    p1 = Tolb*[+Lo -lo 1]';
    p2 = Tolb*[-Lo -lo 1]';
    p3 = Tolb*[+Lo +lo 1]';
    p4 = Tolb*[-Lo +lo 1]';
    line([p1(1) p2(1)],[p1(2) p2(2)],'LineWidth',2,'MarkerSize',10,'color',[0 0 1])
    line([p1(1) p3(1)],[p1(2) p3(2)],'LineWidth',2,'MarkerSize',10,'color',[0 0 1])
    line([p2(1) p4(1)],[p2(2) p4(2)],'LineWidth',2,'MarkerSize',10,'color',[0 0 1])
    line([p3(1) p4(1)],[p3(2) p4(2)],'LineWidth',2,'MarkerSize',10,'color',[0 0 1])
end

function Dibujar_Sistema_Referencia_3D (aTb,s)
    apx = aTb*[0.1 0 0 1]';
    apy = aTb*[0 0.1 0 1]';
    apz = aTb*[0 0 0.1 1]';
    atb = aTb(1:3,4);
    
    line([atb(1) apx(1)],[atb(2) apx(2)],[atb(3) apx(3)],'color',[1 0 0],'LineWidth',3)
    line([atb(1) apy(1)],[atb(2) apy(2)],[atb(3) apy(3)],'color',[0 1 0],'LineWidth',3)
    line([atb(1) apz(1)],[atb(2) apz(2)],[atb(3) apz(3)],'color',[0 0 1],'LineWidth',3)
    
    plot3(atb(1),atb(2),atb(3),'k.','MarkerSize',25)
    text(atb(1),atb(2),atb(3)-0.075,s,'Color','blue','FontSize',15)
end

function z = Penalty (x,xl,xu)
    z = 0;
    m = numel(xl);
    
    for i=1:m
        if xl(i)<x(i)
            z = z + 0;
        else
            z = z + (x(i)-xl(i))^2;
        end

        if x(i)<xu(i)
            z = z + 0;
        else
            z = z + (x(i)-xu(i))^2;
        end
    end
end