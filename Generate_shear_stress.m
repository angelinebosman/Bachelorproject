r1 = 1; %straal buitencirkel
r2 = 0.85; %straal binnencirkel 
n = 200; %hoeveelheid pixels

Shear_stress =sparse(n,n);

for j= 1:n
    for i=1:n
        [placement, edge] = indices(i,j,2*r1/n,r1,r2); %bepaald waar in het domein ligt het roosterpunt
        if placement == "binnenrand"
            if i < 100
                Shear_stress(i,j) = 0.001;
            else
                Shear_stress(i,j) = 1;
            end
        end
    end
end

figure;
hold on;
title('Generated Shear stress');
imagesc(Shear_stress); %display image with scaled colors
colormap('gray') 
axis('square');
axis off;
colorbar;
