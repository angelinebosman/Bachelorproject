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
    S = linspace(0.1,5,20);
    size_lumen = zeros(10,2);
    sum_m = zeros(20);
    
%% Lees de concentratie van LDL op t=0
    fn = double(imread('Beginplot_symmetrisch2.gif')) / 255; fn(fn < 0) = 0;
	K = abs(imresize(fn, [n n])); %zorgt dat het plaatje de juiste afmeting heeft
    b = sparse(reshape(K,[n*n,1]) /1000);
    
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