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

% qa = [0.0 0.5 0.0 0 pi/2 -pi/4 pi/4 0.0]';
% qb = [0.0 -0.5 0.0 0 pi/2 -pi/4 pi/4 0.0]';

qa = [0.3 0.5 pi/4 0 pi/2 0 pi/2 pi/4]';
qb = [0.3 -0.5 -pi/4 0 pi/2 0 pi/2 -pi/4]';

Tx = Txy(qa(1));
Ty = Tx*Tyb(qa(2));
Tb = Ty*Tb0(qa(3));
T1 = Tb*T01(qa(4));
T2 = T1*T12(qa(5));
T3 = T2*T23(qa(6));
T4 = T3*T34(qa(7));
T5 = T4*T45(qa(8));
wTe = T5;

wTe
Dibujar_MM(qa,Tx,Ty,Tb,T1,T2,T3,T4,T5,'a')


Tx = Txy(qb(1));
Ty = Tx*Tyb(qb(2));
Tb = Ty*Tb0(qb(3));
T1 = Tb*T01(qb(4));
T2 = T1*T12(qb(5));
T3 = T2*T23(qb(6));
T4 = T3*T34(qb(7));
T5 = T4*T45(qb(8));
wTe = T5;

wTe
Dibujar_MM(qb,Tx,Ty,Tb,T1,T2,T3,T4,T5,'b')



%% Funciones
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
