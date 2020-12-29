function [x,y,z] = Trayectoria (trac,N)
    switch trac
        case 1
            x = 0.5*ones(1,N);
            y = linspace(0.5,1.5,N);
            z = 0.4 + 0.05*sin(30*y);
        case 2
            theta = linspace(0,2*pi,N);

            x = 0.5*ones(1,N);
            y = 0.0 + 0.05*cos(theta);
            z = 0.4 + 0.05*sin(theta);
        case 3
            x = 0.5*ones(1,N);
            y = linspace(0.5,1.5,N);
            z = 0.4 + 0.1*sin(30*y);
            
            I = z>0.45;
            z(I) =  0.45*ones(1,numel(z(I)));

            I = z<0.35;
            z(I) =  0.35*ones(1,numel(z(I)));
        case 4
            theta = linspace(0,2*pi,N);
            r = 0.035 + 0.015*cos(3*theta);
            
            x = 0.5*ones(1,N);
            y = 0.0 + r.*cos(theta);
            z = 0.4 + r.*sin(theta);
    end