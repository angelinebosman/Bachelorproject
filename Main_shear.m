clear all;
%% initialisatie parameters   
    r1 = 1; %straal buitencirkel
    r2 = 0.85; %straal binnencirkel 
    n = 200; %hoeveelheid pixels
    l = 2*r1/n; %stapgrootte
    dt = 86400*365; %tijdstapgrootte (= 1 jaar)
    k = 10; %aantal jaren
    Z = zeros(n,n);
    O = ones(n,n);
    S = [0.1,1, 2.5, 5];
    size_lumen = zeros(10,4);
    
%% Lees de concentratie van LDL op t=0
    fn = double(imread('Beginplot_symmetrisch2.gif')) / 255; fn(fn < 0) = 0;
	K = abs(imresize(fn, [n n])); %zorgt dat het plaatje de juiste afmeting heeft
    b = sparse(reshape(K,[n*n,1]) /1000);
    
    
%% Vind de oplossingen van de discretisatie systemen
for j=1:4
    shear_stress = S(j)*ones(n,n);
	L = LDL_metnoord(n,r1,r2,Z,Z,b,'n','n');L(isnan(L))=0; 
	C = Chemoattractant_metnoord(n,r1,r2,O,O,O,'n');C(isnan(C))=0;
	m = Monocyten_metnoord(n,r1,r2,C,L,shear_stress,'n');m(isnan(m))=0;
	M = Macrofagen_metnoord(n,r1,r2,m,C,L,'n');M(isnan(M))=0;
	F = Foam_metnoord(n,r1,r2,M,L,'n');F(isnan(F))=0;
    
%% Bereken de grote van het lumen
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
    xlabel('Tijd (in jaren)')
    ylabel('Oppervlakte lumen (in mm)')
    hold on
    plot(tijdstappen,size_lumen(:,2))
    hold on
    plot(tijdstappen,size_lumen(:,3))
    hold on
    plot(tijdstappen,size_lumen(:,4))
    hold off
    



