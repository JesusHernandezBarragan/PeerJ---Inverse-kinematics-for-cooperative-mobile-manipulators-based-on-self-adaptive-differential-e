clear all
close all
clc


N = 200;

[x,y,z] = Trayectoria(1,N);
figure
hold on
grid on
plot3(x,y,z,'*')
axis equal
view([-70 20])

[x,y,z] = Trayectoria(2,N);
figure
hold on
grid on
plot3(x,y,z,'*')
axis equal
view([-70 20])

[x,y,z] = Trayectoria(3,N);
figure
hold on
grid on
plot3(x,y,z,'*')
axis equal
view([-70 20])

[x,y,z] = Trayectoria(4,N);
figure
hold on
grid on
plot3(x,y,z,'*')
axis equal
view([-70 20])
