clear all;
% Menu
   title = 'Which version of the model do you want to run? (1 for yes 0 for no)';
   prompt = { 'Differentiating shear stress', 'moving boundary problem'};
   defaultanswer={'0','0'};
   params=inputdlg(prompt,title, [1 110], defaultanswer);
   shear            = str2num(char(params(1)));
   moving_boundary  = str2num(char(params(2)));

r1 = 1; %straal buitencirkel
r2 = 0.85; %straal binnencirkel 
n = 200; %hoeveelheid pixels
l = 2*r1/n; %stapgrootte
dt = 0.5; %tijdstapgrootte

%do you want to have differentiating shear stress:
%shear = 1; %1 for yes, 0 for no

%do you want to solve the moving boundary problem:
%moving_boundary = 0; %1 for yes, 0 for no


%generate the shear stress over the domain
shear_stress =ones(n,n);
if shear == 1
    for j= 1:n
        for i=1:n
            [placement, edge] = indices(i,j,l,r1,r2); %bepaald waar in het domein ligt het roosterpunt
            if placement == "binnenrand"
                if i < 100
                    if j < 100
                        shear_stress(i,j) = 0.01;
                    end
                end
            end
        end
    end
end


    fn = double(imread('Beginplot_plaque_links2.gif')) / 255;
    K = imresize(fn, [n n]); %zorgt dat het plaatje de juiste afmeting heeft
    b = 5*reshape(K,[n*n,1]);
    
    Z = zeros(n,n);
    L = LDL(n,r1,r2,Z,Z,b);L(isnan(L))=0;
    C = Chemoattractant(n,r1,r2,Z,Z,Z);C(isnan(C))=0;
    
    m = Monocyten(n,r1,r2,C,L,shear_stress,'n');m(isnan(m))=0;
    M = Macrofagen(n,r1,r2,m,C,L);M(isnan(M))=0;
    F = Foam(n,r1,r2,M,L);F(isnan(F))=0;
    L = LDL(n,r1,r2,m,M,b);L(isnan(L))=0;
    C = Chemoattractant(n,r1,r2,m,M,F);C(isnan(C))=0;
    
    if moving_boundary == 1
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
        Levelset = LevelSet(n,r1,r2,dt,L,m,M,p);
        for i=1:20
            Levelset = LevelSet(n,r1,r2,dt,L,m,M,Levelset);
        end
        [xv,yv] = get_xv_and_yv(Levelset);
        %m = Monocyten(n,r1,r2,C,L,shear_stress,'y',xv,yv);
    end

%plaatjes maken
figure;hold on; title('Monocyten');imagesc(m);colormap('gray');axis off;colorbar;
figure;hold on;title('Macrofagen');imagesc(M);colormap('gray') ;axis off;colorbar;
figure;hold on;title('Chemoattractant');imagesc(C); colormap('gray');axis off;colorbar;
figure;hold on;title('LDL');imagesc(L);colormap('gray');axis off;colorbar;
figure;hold on;title('Foam cells');imagesc(F);colormap('gray');axis off;colorbar;
%figure;hold on;title('Nieuwe binnenrand');imagesc(Levelset);colormap('gray');axis off;colorbar;











