clc; clear all ;
M = 9;
N = 100;
r1 = 0.85 ; % inner radius 
r2 = 1 ;  % outer radius
[T,R] = meshgrid(linspace(0,2*pi,N),linspace(r1,r2,M));
% Convert grid to cartesian coordintes
X = R.*cos(T); 
Y = R.*sin(T);
T = delaunay(X,Y);
figure;
hold on;
trimesh(T,X,Y);

%patch middle
t = linspace(0, 2*pi);
X1 = r1*cos(t);
Y1 = r1*sin(t);
patch(X1,Y1,'white','EdgeColor','none')
set(gca,'XColor', 'none','YColor','none') %removes x and y axis
hold off
