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

q0 = [0.0 0.5 pi/4 0 pi/2 -pi/4 pi/4 0.0 0.0 -0.5 pi/4 0 pi/2 -pi/4 pi/4 0.0]'; 

td1 = [0.7 0.0 0.5]';
prd = [0.0 -1.0 0.0]';
td2 = td1 + prd;

N = 200;
[x,y,z] = Trayectoria(4,N);
Q = zeros(16,N);
time = zeros(1,N);
fitness = zeros(1000,N); % 1000

for i=1:N
    tic
    
    td1 = [x(i) y(i) z(i)]';
    td2 = td1 + prd;
    
    [xgb,fitness(:,i)] = IK (q0,td1,td2,Txy,Tyb,Tb0,T01,T12,T23,T34,T45);
    q0 = xgb;
    
    Q(:,i) = xgb;
    
    time(i) = toc;
    disp(['iter = ' num2str(i) ', time = ' num2str(time(i))])
end

save('CMAES_4','x','y','z','prd','Q','N','time','fitness')
return

figure
hold on
grid on

for i=1:N
    clf
    
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

    Dibujar_MM(q1,Tx,Ty,Tb,T1,T2,T3,T4,T5,'a')

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

    Dibujar_MM(q2,Tx,Ty,Tb,T1,T2,T3,T4,T5,'b')

%     plot3(td1(1),td1(2),td1(3),'go','MarkerSize',15,'LineWidth',2)
%     plot3(td1(1),td1(2),td1(3),'rx','MarkerSize',15,'LineWidth',2)
% 
%     plot3(td2(1),td2(2),td2(3),'go','MarkerSize',15,'LineWidth',2)
%     plot3(td2(1),td2(2),td2(3),'rx','MarkerSize',15,'LineWidth',2)

    view([-21 50])
    
    pause(0.01)
end


%% Funciones
function [xgb,plot_fit] = IK (q0,td1,td2,Txy,Tyb,Tb0,T01,T12,T23,T34,T45)
    dim = 16;
    gen = 1000; % 1000
    tol = 0.0001;
    lambda = 30; %30
    
    lb = [-1.5 -1.5 -pi -169*(pi/180) -65*(pi/180) -150*(pi/180) -102.2*(pi/180) -167.5*(pi/180) -1.5 -1.5 -pi -169*(pi/180) -65*(pi/180) -150*(pi/180) -102.2*(pi/180) -167.5*(pi/180)]';
    ub = [+1.5 +1.5 +pi +169*(pi/180) +90*(pi/180) +146*(pi/180) +102.5*(pi/180) +167.5*(pi/180) +1.5 +1.5 +pi +169*(pi/180) +90*(pi/180) +146*(pi/180) +102.5*(pi/180) +167.5*(pi/180)]';

    % Selection and Recombination:
    % lambda = round(4+(3*log(dim)));
    mu = round(lambda/2);
    w = log(mu+1/2)-log(1:mu)';
    w = w/sum(w);
    mueff = 1/sum(w.^2);

    % Step-size control:
    cs = (mueff+2)/(dim+mueff+5);
    ds = 1+2*max(0,sqrt((mueff-1)/(dim+1))-1)+cs;

    % Covariance matrix adaptation:
    cc = (4+mueff/dim)/(dim+4+2*mueff/dim);
    c1 = 2/((dim+1.3)^2+mueff);
    cmu = min(1-c1,2*(mueff-2+1/mueff)/((dim+2)^2+mueff));
    E = (dim^0.5)*(1-1/(4*dim)+1/(21*dim^2));

    %% Initialization
    ps = zeros(dim,1);
    pc = zeros(dim,1);

    B = eye(dim);
    D = eye(dim);
    C = B*D*B';
    d = D;

    for i = 1:dim
        d(i,i) = d(i,i)^-1;
    end

    invsqrtC = B*d*B';

    xmean = rand(dim,1);
    sigma = 0.5;

    %% Until Termination Criterion Met
    % CMA-ES Algorithm
    fit = zeros(lambda,1);
    xk = zeros(dim,lambda);
    yk = zeros(dim,lambda);
    
    f = @(t1,t2,q) 1.5*(norm(td1-t1)+norm(td2-t2)) + 0.6*norm(q0-q) + 100*Penalty(q,lb,ub);
	
    plot_fit = zeros(1,gen);
    
    for g = 1:gen
        % Sample new population of search points, for k = 1,...,lambda
        for k = 1:lambda
            zk = randn(dim,1);
            yk(:,k) = B*D*zk;
            xk(:,k) = xmean+sigma*yk(:,k);

            q1 = xk(1:8,k);
            Tx = Txy(q1(1));
            Ty = Tx*Tyb(q1(2));
            Tb = Ty*Tb0(q1(3));
            T1 = Tb*T01(q1(4));
            T2 = T1*T12(q1(5));
            T3 = T2*T23(q1(6));
            T4 = T3*T34(q1(7));
            T5 = T4*T45(q1(8));
            Tq1 = T5;

            q2 = xk(9:16,k);
            Tx = Txy(q2(1));
            Ty = Tx*Tyb(q2(2));
            Tb = Ty*Tb0(q2(3));
            T1 = Tb*T01(q2(4));
            T2 = T1*T12(q2(5));
            T3 = T2*T23(q2(6));
            T4 = T3*T34(q2(7));
            T5 = T4*T45(q2(8));
            Tq2 = T5;

            fit(k) = f(Tq1(1:3,4),Tq2(1:3,4),xk(:,k)); 
        end

        % Selection and recombination
        [fit,indfit] = sort(fit);
        xold = xmean;

        xmean = xk(:,indfit(1:mu))*w;
        yw = (xmean-xold)/sigma;

        % Step-size control 
        ps = (1-cs)*ps+sqrt(cs*(2-cs)*mueff)*invsqrtC*yw;
        sigma = sigma*exp((cs/ds)*(norm(ps)/E-1));

        % Covariance matrix adaptation
        if (norm(ps)/sqrt(1-(1-cs)^(2*(g+1)))) < ((1.4+2/(dim+1))*E)
            hs = 1;
        else
            hs = 0;
        end

        pc = (1-cc)*pc+hs*sqrt(cc*(2-cc)*mueff)*yw;
        dhs = (1-hs)*cs*(2-cc);
        C = (1-c1-cmu)*C+c1*(pc*pc'+dhs*C)+cmu*(yk(:,indfit(1:mu))*diag(w)*yk(:,indfit(1:mu))');

        C = triu(C)+triu(C,1)';
        [B,D] = eig(C);
        B = real(B);
        D = real(D);  
        D = real(diag(sqrt(diag(D))));

        for i = 1:dim
            d(i,i) = d(i,i)^-1;
        end

        invsqrtC = B*d*B';

        q1 = xk(1:8,indfit(1));
        Tx = Txy(q1(1));
        Ty = Tx*Tyb(q1(2));
        Tb = Ty*Tb0(q1(3));
        T1 = Tb*T01(q1(4));
        T2 = T1*T12(q1(5));
        T3 = T2*T23(q1(6));
        T4 = T3*T34(q1(7));
        T5 = T4*T45(q1(8));
        Tq1 = T5;

        q2 = xk(9:16,indfit(1));
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

       % plot_fit(g) = fit(1);
        plot_fit(g) = e1+e2;
        
        if(e1<tol && e2<tol)
            break
        end
    end

    xgb = xk(:,indfit(1));
%     fmin = fit(1);
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



