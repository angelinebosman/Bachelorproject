%% assign parameter values
r1 = 1; %straal buitencirkel
r2 = 0.85; %straal binnencirkel 
n = 200; %hoeveelheid pixels
l = 2*r1/n; %stapgrootte
dt = 86400*30*12; %tijdstapgrootte (= 12 maanden)
k = 20; %aantal jaren
L(L < 0) = 0;
m(m < 0) = 0;
M(M < 0) = 0;

%% create the approximated initial binnenrand equation
 phi_0 = zeros(n,n); 
for i=1:n
     for j=1:n
         [placement, edge] = indices(i,j,2*r1/n,r1,r2);
         if placement == "binnenrand"
             phi_0(i,j) = 1;
         end
     end
end

%% apply level set method
phi = LevelSet_constante_S(n,r1,dt,S,phi);phi(phi < 0) = 0; 
for i=1:k-1
    phi(phi > 0) = 1;phi = LevelSet_constante_S(n,r1,dt,L,m,M,phi); phi(phi < 0) = 0; 
end
%% create plots
[xv,yv] = get_xv_and_yv(phi);
figure
plot(xv,yv,':r','LineWidth',2.1)
hold on
set(gca,'XColor', 'none','YColor','none') %removes x and y axis
[xc,yc] = get_xv_and_yv(phi_0);
plot(xc,yc,':b','LineWidth',2)
text = sprintf('Oplossing level set methode voor constante S na %0.1f jaar', (k+1)/2);
title(text)
legend('nieuwe binnenrand', 'oude binnenrand') 