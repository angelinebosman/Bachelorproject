%% assign parameter values
r1 = 1; %straal buitencirkel
r2 = 0.85; %straal binnencirkel 
n = 200; %hoeveelheid pixels
l = 2*r1/n; %stapgrootte
dt = 86400*30; %tijdstapgrootte (86300 sec = 1 dag)
k = 5; %aantal maanden
O = ones(n,n);

%% create the approximated initial binnenrand equation
 phi_0 = zeros(n,n); 
for i=1:n
     for j=1:n
         [placement, edge] = indices(i,j,2*r1/n,r1,r2);
         if placement == "binnenrand"
             phi_0(i,j) = 0.1;
         end
     end
end

%% apply level set method
Levelset = LevelSet(n,r1,dt,O,O,O,phi_0);
Levelset(Levelset < 0) = 0; %zit geen rand
for i=1:k
    Levelset = LevelSet(n,r1,dt,O,O,O,Levelset);
    Levelset(Levelset < 0) = 0; %zit geen rand
end

%% create plots
[xv,yv] = get_xv_and_yv(Levelset);
figure
plot(xv,yv,'r+')
hold on
%set(gca,'XColor', 'none','YColor','none') %removes x and y axis
[xc,yc] = get_xv_and_yv(phi_0);
plot(xc,yc,'bo')
legend('nieuwe binnenrand', 'oude binnenrand') 