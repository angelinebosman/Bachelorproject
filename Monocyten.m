function [output] = Monocyten(n,r1,r2,C,L,shear_stress,moving_bdry,xv, yv)
    if (~exist('xv', 'var')) %createert een optionele input parameter
        opt = true;
    end
    if (~exist('yv', 'var'))
        opt = true;
    end
    %% Parameters
    mu_1 = 0.05;
    d1 = 0.1; %diffusieconstante bloedwand naar lumen
    xi_1 = 10; %chemotactische gevoeligheid
    rho_0 = 5*10^(-7);%1*10^(-4);
    rho_ster = 1;
    rho_m = 125;%0.125;%7;  %Omzetting van monocyten naar macrofagen

    %% Maak A en f
    l = 2*r1/n; %lengte tussen roosterpunten: diameter domein gedeeld door aantal roosterpunten
    A = sparse(n^2,n^2); %matrix die alleen de niet-nul waardes opslaat
    f = sparse(n^2,1);
    for j= 1:n
        for i=1:n
            k = i+(j-1)*n; %horizontale nummering
            if moving_bdry == 'y'
                [placement, edge] = indices_moving_bdry(i,j,l, r1, xv, yv);
            else
                [placement, edge] = indices(i,j,l,r1,r2); %bepaald waar in het domein het roosterpunt ligt
            end
            if placement == "inside"
                A(k,k) = -4*mu_1 + xi_1/C(i,j) * rho_0/rho_ster * (4*C(i,j) - C(i-1,j) - C(i,j-1) - C(i,j+1) - C(i+1,j)) - l^2/4*rho_m;
                A(k+n,k) = mu_1; %oost roosterpunt
                A(k-n,k) = mu_1; %west roosterpunt
                A(k-1,k) = mu_1; %zuid roosterpunt
                A(k+1,k) = mu_1; %noord roosterpunt
            elseif placement == "buitenrand"
                if edge == "west"
                    A(k,k) = -3*mu_1 + (xi_1/C(i,j))*(rho_0/rho_ster)*(3*C(i,j) - C(i-1,j) - C(i,j+1)- C(i+1,j)) - l^2/4*rho_m;
                    A(k+n,k) = mu_1; %oost roosterpunt
                    A(k-1,k) = mu_1; %zuid roosterpunt
                    A(k+1,k) = mu_1; %noord roosterpunt
                elseif edge == "oost"
                    A(k,k) = -3*mu_1 + xi_1/C(i,j) * rho_0/rho_ster * (3*C(i,j) - C(i,j-1) - C(i-1,j)- C(i+1,j))   - l^2/4*rho_m;
                    A(k-n,k) = mu_1; %west roosterpunt
                    A(k-1,k) = mu_1; %zuid roosterpunt
                    A(k+1,k) = mu_1; %noord roosterpunt
                elseif edge == "noord"
                    A(k,k) = -3*mu_1 + xi_1/C(i,j) * rho_0/rho_ster * (3*C(i,j) - C(i,j-1) - C(i,j+1) - C(i-1,j))   - l^2/4*rho_m;
                    A(k-n,k) = mu_1; %west roosterpunt
                    A(k+n,k) = mu_1; %oost roosterpunt
                    A(k-1,k) = mu_1; %zuid roosterpunt
                elseif edge == "zuid"
                    A(k,k) = -3*mu_1 + xi_1/C(i,j) * rho_0/rho_ster * (3*C(i,j) - C(i,j+1) - C(i,j-1) - C(i+1,j))- l^2/4*rho_m;
                    A(k-n,k) = mu_1; %west roosterpunt
                    A(k+n,k) = mu_1; %oost roosterpunt
                    A(k+1,k) = mu_1; %noord roosterpunt
                elseif edge == "oost_zuid"
                    A(k,k) = -2*mu_1 + xi_1/C(i,j) * rho_0/rho_ster * (2*C(i,j) - C(i,j-1) - C(i+1,j)) - l^2/4*rho_m;
                    A(k-n,k) = mu_1; %west roosterpunt
                    A(k+1,k) = mu_1; %noord roosterpunt
                elseif edge == "west_zuid"
                    A(k,k) = -2*mu_1 + xi_1/C(i,j) * rho_0/rho_ster * (2*C(i,j) - C(i,j+1) - C(i+1,j)) - l^2/4*rho_m;
                    A(k+n,k) = mu_1; %oost roosterpunt
                    A(k+1,k) = mu_1; %noord roosterpunt
                elseif edge == "west_noord"
                    A(k,k) =  -2*mu_1 + xi_1/C(i,j) * rho_0/rho_ster * (2*C(i,j) - C(i,j+1) - C(i-1,j))- l^2/4*rho_m;
                    A(k+n,k) = mu_1; %oost roosterpunt
                    A(k-1,k) = mu_1; %zuid roosterpunt
                elseif edge == "oost_noord"
                    A(k,k) =  -2*mu_1 + xi_1/C(i,j) * rho_0/rho_ster * (2*C(i,j) - C(i-1,j) - C(i,j-1)) - l^2/4*rho_m;
                    A(k-n,k) = mu_1; %west roosterpunt
                    A(k-1,k) = mu_1; %zuid roosterpunt
                end
            elseif placement == "binnenrand" 
                f(k,1) = -mu_1*(l^2/4)* Lambda(shear_stress(i,j), L(i,j));
                if edge == "west"
                    A(k,k) = -3*mu_1 + xi_1/C(i,j) * rho_0/rho_ster * (3*C(i,j) - C(i,j-1) - C(i-1,j)  - C(i+1,j))   - l^2/4*rho_m - sqrt(3)/3*mu_1*d1;
                    A(k-n,k) = mu_1; %west roosterpunt
                    A(k-1,k) = mu_1; %zuid roosterpunt
                    A(k+1,k) = mu_1; %noord roosterpunt
                elseif edge == "oost"
                    A(k,k) = -3*mu_1 + xi_1/C(i,j) * rho_0/rho_ster * (3*C(i,j) - C(i-1,j) - C(i,j+1)  - C(i+1,j))   - l^2/4*rho_m  - sqrt(3)/3*mu_1*d1;
                    A(k+n,k) = mu_1; %oost roosterpunt
                    A(k-1,k) = mu_1; %zuid roosterpunt
                    A(k+1,k) = mu_1; %noord roosterpunt
                elseif edge == "noord"
                    A(k,k) = -3*mu_1 + xi_1/C(i,j) * rho_0/rho_ster * (3*C(i,j) - C(i,j+1) - C(i,j-1) - C(i+1,j))   - l^2/4*rho_m  - sqrt(3)/3*mu_1*d1;
                    A(k+n,k) = mu_1; %oost roosterpunt
                    A(k-n,k) = mu_1; %west roosterpunt
                    A(k+1,k) = mu_1; %noord roosterpunt
                elseif edge == "zuid"
                    A(k,k) = -3*mu_1 + xi_1/C(i,j) * rho_0/rho_ster * (3*C(i,j) - C(i,j-1) - C(i,j+1) - C(i-1,j))   - l^2/4*rho_m - sqrt(3)/3*mu_1*d1;
                    A(k+n,k) = mu_1; %oost roosterpunt
                    A(k-n,k) = mu_1; %west roosterpunt
                    A(k-1,k) = mu_1; %zuid roosterpunt
                elseif edge == "west_noord"
                    A(k,k) = -2*mu_1 + xi_1/C(i,j) * rho_0/rho_ster * (2*C(i,j) - C(i,j-1) - C(i+1,j)) - l^2/4*rho_m - sqrt(3)/3*mu_1*d1;
                    A(k-n,k) = mu_1; %west roosterpunt
                    A(k+1,k) = mu_1; %noord roosterpunt
                elseif edge == "oost_noord"
                    A(k,k) = -2*mu_1 + xi_1/C(i,j) * rho_0/rho_ster * (2*C(i,j) - C(i+1,j) - C(i,j+1)) - l^2/4*rho_m - sqrt(3)/3*mu_1*d1;
                    A(k+n,k) = mu_1; %oost roosterpunt
                    A(k+1,k) = mu_1; %noord roosterpunt
                elseif edge == "west_zuid"
                    A(k,k) = -2*mu_1 + xi_1/C(i,j) * rho_0/rho_ster * (2*C(i,j) - C(i,j-1) - C(i-1,j))   - l^2/4*rho_m - sqrt(3)/3*mu_1*d1;
                    A(k-n,k) = mu_1; %west roosterpunt
                    A(k-1,k) = mu_1; %zuid roosterpunt
                elseif edge == "oost_zuid"
                    A(k,k) = -2*mu_1 + xi_1/C(i,j) * rho_0/rho_ster * (2*C(i,j) - C(i-1,j) - C(i,j+1))   - l^2/4*rho_m - sqrt(3)/3*mu_1*d1;
                    A(k+n,k) = mu_1; %oost roosterpunt
                    A(k-1,k) = mu_1; %zuid roosterpunt
                end
            end
        end
    end

    %% Vind de oplossing voor Au=f
    u = A\f;u(isinf(u))=0;
    output = reshape(u,[n n])
end