r1 = 1; %straal buitencirkel
r2 = 0.85; %straal binnencirkel 
n = 200; %hoeveelheid pixels
l = 2*r1/n; %stapgrootte
dt = 5; %tijdstapgrootte
O = ones(n,n);

 p = sparse(n,n); %create the plot of the initial binnenrand equation
 for i=1:n
     for j=1:n
         [placement, edge] = indices(i,j,2*r1/n,r1,r2);
         if placement == "binnenrand"
             p(i,j) = 1;
         end
     end
 end
p = reshape(p,[n*n,1]);
Levelset = LevelSet(n,r1,r2,dt,O,O,O,p);
Levelset(isnan(Levelset))=0;
for i=1:5
    Levelset = LevelSet(n,r1,r2,dt,O,O,O,Levelset);
    Levelset(isnan(Levelset))=0;
end
imagesc(Levelset)