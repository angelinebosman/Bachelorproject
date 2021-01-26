clear all;

%% Menu
   Title = 'Welke versie van het model wil je runnen? (1 voor ja 0 voor nee)';
   prompt = { 'Plaatsafhankelijke shear stress', 'Bewegend-rand probleem', 'Homogeen probleem','Monocyt dichtheid bij wisselende shear stress', 'Oppervlakte lumen bij wisselende shear stress'};
   defaultanswer={'0','0','0','0', '0'};
   params=inputdlg(prompt,Title, [1 110], defaultanswer);
   shear            = str2num(char(params(1)));
   bewegende_rand   = str2num(char(params(2)));
   homogeen         = str2num(char(params(3)));
   optie1           = str2num(char(params(4)));
   optie2           = str2num(char(params(5)));
   
   
%% initialisatie parameters   
    r1 = 1; %straal buitencirkel
    r2 = 0.85; %straal binnencirkel 
    n = 200; %hoeveelheid pixels
    l = 2*r1/n; %stapgrootte
    dt = 86400*365; %tijdstapgrootte (= 1 jaar)
    k = 10; %aantal jaren
    Z = zeros(n,n);
    O = ones(n,n);
    
    
%% Lees de concentratie van LDL op t=0
	fn = double(imread('Beginplot_symmetrisch2.gif')) / 255; fn(fn < 0) = 0;
	K = abs(imresize(fn, [n n])); %zorgt dat het plaatje de juiste afmeting heeft
	b = sparse(reshape(K,[n*n,1]) /1000);
    
        
%% Maak een plot van de totale monocyt dichtheid bij wisselende shear stress
    if optie1 == 1
        size_lumen = zeros(10,2);
        sum_m = zeros(20);
        L = LDL_metnoord(n,r1,r2,Z,Z,b,'n','n');L(isnan(L))=0; 
        C = Chemoattractant_metnoord(n,r1,r2,O,O,O,'n');C(isnan(C))=0;
        for j=1:20
            shear_stress = S(j)*ones(n,n);
            m = Monocyten_metnoord(n,r1,r2,C,L,shear_stress,'n');m(isnan(m))=0;
            sum_m(j) = sum(sum(m));
        end
        plot(S,sum_m)
        xlabel("Shear stress parameter \sigma_w")
        ylabel("Som van de monocyten dichtheid in de bloedwand")
        return %run niet meer de rest van het programma
    end

    
 %% Maak een plot van het oppervlakte van het lumen bij wisselende shear stress
    if optie2 == 1
        S = [0.1,1, 2.5, 5];
        size_lumen = zeros(10,4);
        for j=1:4
            shear_stress = S(j)*ones(n,n);
            L = LDL_metnoord(n,r1,r2,Z,Z,b,'n','n');L(isnan(L))=0; 
            C = Chemoattractant_metnoord(n,r1,r2,O,O,O,'n');C(isnan(C))=0;
            m = Monocyten_metnoord(n,r1,r2,C,L,shear_stress,'n');m(isnan(m))=0;
            M = Macrofagen_metnoord(n,r1,r2,m,C,L,'n');M(isnan(M))=0;
            F = Foam_metnoord(n,r1,r2,M,L,'n');F(isnan(F))=0;
            phi = zeros(n,n); %phi_t is phi op tijdstap t
            for l=1:n
                 for p=1:n
                    [placement, edge] = indices(l,p,2*r1/n,r1,r2);
                    if placement == "binnenrand"
                    phi(l,p) = 1;
                    end
                 end
            end
            t = 0;
            [xv,yv] = get_xv_and_yv(phi);
            size_lumen(1,j) = polyarea(xv,yv);
            phi = LevelSet(n,r1,dt,L,m,M,phi);phi(phi < 0) = 0; %phi op tijdstap t+1
            for i=2:k
                [xv,yv] = get_xv_and_yv(phi);
                size_lumen(i,j) = polyarea(xv,yv);
                if mod(k,5) == 0
                    L = LDL_metnoord(n,r1,r2,m,M,b,'n','y',xv,yv);L(isnan(L))=0;
                    C = Chemoattractant_metnoord(n,r1,r2,m,M,F,'y',xv,yv);C(isnan(C))=0;
                    m = Monocyten_metnoord(n,r1,r2,C,L,shear_stress,'y',xv,yv);m(isnan(m))=0;
                    M = Macrofagen_metnoord(n,r1,r2,m,C,L,'y',xv,yv);M(isnan(M))=0;
                    F = Foam_metnoord(n,r1,r2,M,L,'y',xv,yv);F(isnan(F))=0;
                end
                phi = LevelSet(n,r1,dt,L,m,M,phi); phi(phi < 0) = 0;
            end         
        end
        tijdstappen = 1:k;
        plot(tijdstappen,size_lumen(:,1))
        xlabel('Tijd (in jaren)');ylabel('Oppervlakte lumen (in mm)');
        hold on
        plot(tijdstappen,size_lumen(:,2));hold on
        plot(tijdstappen,size_lumen(:,3));hold on
        plot(tijdstappen,size_lumen(:,4));
        legend('\sigma_w = 0.1','\sigma_w = 1','\sigma_w = 2.5', '\sigma_w = 5')
        hold off
        return %run niet meer de rest van het programma
    end

    
%% genereer de plaatsafhankelijke shear stress
    shear_stress = ones(n,n);
    if shear == 1
        for j= 1:n
            for i=1:n
                    if i < 100
                        if j < 100
                            shear_stress(i,j) = 0.001;
                        end
                    end
            end
        end    
    end
    
    
%% Homogeen probleem
    if homogeen == 1
        L = LDL_metnoord(n,r1,r2,O,O,Z,'y','n');L(isnan(L))=0; 
        C = Chemoattractant_metnoord(n,r1,r2,O,O,O,'n');C(isnan(C))=0;
        m = Monocyten_metnoord(n,r1,r2,C,L,shear_stress,'n');m(isnan(m))=0;
        M = Macrofagen_metnoord(n,r1,r2,C,m,L,'n');M(isnan(M))=0;
        F = Foam_metnoord(n,r1,r2,M,L,'n');F(isnan(F))=0;
        C = Chemoattractant_metnoord(n,r1,r2,m,M,F,'n');C(isnan(C))=0;
    else         
%% Vind de oplossingen van de discretisatie systemen    
        L = LDL_metnoord(n,r1,r2,Z,Z,b,'n','n');L(isnan(L))=0; 
        C = Chemoattractant_metnoord(n,r1,r2,O,O,O,'n');C(isnan(C))=0;
        m = Monocyten_metnoord(n,r1,r2,C,L,shear_stress,'n');m(isnan(m))=0;
        M = Macrofagen_metnoord(n,r1,r2,m,C,L,'n');M(isnan(M))=0;
        F = Foam_metnoord(n,r1,r2,M,L,'n');F(isnan(F))=0;
        C = Chemoattractant_metnoord(n,r1,r2,m,M,F,'n');C(isnan(C))=0;
        L = LDL_metnoord(n,r1,r2,m,M,b,'n','n');L(isnan(L))=0; 
   
 %%  Los het bewegengd-randprobleem op
        if bewegende_rand == 1
             phi_t = zeros(n,n); %phi_t is phi op tijdstap t
             for i=1:n
                 for j=1:n
                    [placement, edge] = indices(i,j,2*r1/n,r1,r2);
                    if placement == "binnenrand"
                    phi_t(i,j) = 1;
                    end
                 end
             end

            phi = LevelSet(n,r1,dt,L,m,M,phi_t);phi(phi < 0) = 0; %phi op tijdstap t+1
            for i=1:k-1
                [xv,yv] = get_xv_and_yv(phi);
                L = LDL_metnoord(n,r1,r2,m,M,b,'n','y',xv,yv);L(isnan(L))=0;
                C = Chemoattractant_metnoord(n,r1,r2,m,M,F,'y',xv,yv);C(isnan(C))=0;
                m = Monocyten_metnoord(n,r1,r2,C,L,shear_stress,'y',xv,yv);m(isnan(m))=0;
                M = Macrofagen_metnoord(n,r1,r2,m,C,L,'y',xv,yv);M(isnan(M))=0;
                F = Foam_metnoord(n,r1,r2,M,L,'y',xv,yv);F(isnan(F))=0;
                phi = LevelSet(n,r1,dt,L,m,M,phi); phi(phi < 0) = 0;
            end 
        end
    end

    
%% Maak de plaatjes
    figure;hold on; title('Monocyten');imagesc(m);axis off;colorbar;
    figure;hold on;title('Macrofagen');imagesc(M);axis off;colorbar;
    figure;hold on;title('Chemoattractant');imagesc(C);axis off;colorbar; 
    figure;hold on;title('LDL');imagesc(L);axis off;colorbar;
    figure;hold on;title('Foam cells');imagesc(F);axis off;colorbar;
    
    if bewegende_rand == 1
        [xv,yv] = get_xv_and_yv(phi);
        figure
        plot(yv,xv,':r','LineWidth',2.1)
        hold on
        set(gca,'XColor', 'none','YColor','none') %removes x and y axis
        [xc,yc] = get_xv_and_yv(phi_t);
        plot(xc,yc,':b','LineWidth',2)
        text = sprintf('Oplossing level set methode na %1.0f jaar', k);
        title(text)
        legend('nieuwe binnenrand', 'oude binnenrand')
    end

