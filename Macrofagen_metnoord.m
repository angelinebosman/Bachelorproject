function [output] = Macrofagen_metnoord(n,r1,r2,m,C,L,moving_bdry,xv, yv)
    if (~exist('xv', 'var')) %creates an optional input parameter
        opt = true;
    end
    if (~exist('yv', 'var'))
        opt = true;
    end
    
    %% Parameters
    mu_2 = 0.05;
    xi_2 = 10;
    rho_0 = 5*10^(-7);
    rho_ster = 1;
    rho_m = 125;%0.125; %7 %omzetting van monocyten naar macrofagen
    rho_in = 10^(-2); 
    L_th = 1;
    

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
            if placement ~= "outside"
                f(k,1) = (l^2/4 * rho_m * m(i,j));
            end
            if placement == "inside"
                A(k,k) = -4*mu_2 + (xi_2/C(i,j))*(rho_0/rho_ster)*(4*C(i,j) - C(i-1,j) - C(i,j-1) - C(i,j+1) - C(i+1,j)) - ((l^2/4)*(rho_in * L(i,j))/(1+L(i,j)/L_th));
                A(k+n,k) = mu_2; %oost roosterpunt
                A(k-n,k) = mu_2; %west roosterpunt
                A(k-1,k) = mu_2; %zuid roosterpunt
                A(k+1,k) = mu_2; %noord roosterpunt
            elseif placement == "buitenrand"
                if edge == "west"
                    A(k,k) = -3*mu_2 + (xi_2/C(i,j))*(rho_0/rho_ster)*(3*C(i,j) - C(i-1,j) - C(i,j+1)- C(i+1,j)) - ((l^2/4)*(rho_in * L(i,j))/(1+L(i,j)/L_th));
                    A(k+n,k) = mu_2; %oost roosterpunt
                    A(k-1,k) = mu_2; %zuid roosterpunt
                    A(k+1,k) = mu_2; %noord roosterpunt
                elseif edge == "oost"
                    A(k,k) = -3*mu_2 + (xi_2/C(i,j))*(rho_0/rho_ster)*(3*C(i,j) - C(i,j-1) - C(i-1,j)- C(i+1,j))- ((l^2/4)*(rho_in * L(i,j))/(1+L(i,j)/L_th));
                    A(k-n,k) = mu_2; %west roosterpunt
                    A(k-1,k) = mu_2; %zuid roosterpunt
                    A(k+1,k) = mu_2; %noord roosterpunt
                elseif edge == "noord"
                    A(k,k) = -3*mu_2 + (xi_2/C(i,j)) * rho_0/rho_ster * (3*C(i,j) - C(i,j-1) - C(i,j+1) - C(i-1,j)) -((l^2/4)*(rho_in * L(i,j))/(1+L(i,j)/L_th));
                    A(k-n,k) = mu_2; %west roosterpunt
                    A(k+n,k) = mu_2; %oost roosterpunt
                    A(k-1,k) = mu_2; %zuid roosterpunt
                elseif edge == "zuid" 
                    A(k,k) = -3*mu_2 + (xi_2/C(i,j)) * rho_0/rho_ster * (3*C(i,j) - C(i,j+1) - C(i,j-1) - C(i+1,j))- ((l^2/4)*(rho_in * L(i,j))/(1+L(i,j)/L_th));
                    A(k-n,k) = mu_2; %west roosterpunt
                    A(k+n,k) = mu_2; %oost roosterpunt
                    A(k+1,k) = mu_2; %noord roosterpunt
                elseif edge == "oost_zuid"
                    A(k,k) = -2*mu_2 + (xi_2/C(i,j)) * rho_0/rho_ster * (2*C(i,j) - C(i,j-1) - C(i+1,j)) - ((l^2/4)*(rho_in * L(i,j))/(1+L(i,j)/L_th));
                    A(k-n,k) = mu_2; %west roosterpunt
                    A(k+1,k) = mu_2; %noord roosterpunt
                elseif edge == "west_zuid"
                    A(k,k) = -2*mu_2 + (xi_2/C(i,j)) * rho_0/rho_ster * (2*C(i,j) - C(i,j+1) - C(i+1,j)) - ((l^2/4)*(rho_in * L(i,j))/(1+L(i,j)/L_th));
                    A(k+n,k) = mu_2; %oost roosterpunt
                    A(k+1,k) = mu_2; %noord roosterpunt
                elseif edge == "west_noord"
                    A(k,k) =  -2*mu_2 + (xi_2/C(i,j)) * rho_0/rho_ster * (2*C(i,j) - C(i,j+1) - C(i-1,j)) - ((l^2/4)*(rho_in * L(i,j))/(1+L(i,j)/L_th));
                    A(k+n,k) = mu_2; %oost roosterpunt
                    A(k-1,k) = mu_2; %zuid roosterpunt
                elseif edge == "oost_noord"
                    A(k,k) =  -2*mu_2 + (xi_2/C(i,j)) * rho_0/rho_ster * (2*C(i,j) - C(i-1,j) - C(i,j-1)) - ((l^2/4)*(rho_in * L(i,j))/(1+L(i,j)/L_th));
                    A(k-n,k) = mu_2; %west roosterpunt
                    A(k-1,k) = mu_2; %zuid roosterpunt
                end
            elseif placement == "binnenrand"
                if edge == "west"
                    A(k,k) = -3*mu_2 + (xi_2/C(i,j)) * rho_0/rho_ster * (3*C(i,j) - C(i,j-1) - C(i-1,j)  - C(i+1,j)) - ((l^2/4)*(rho_in * L(i,j))/(1+L(i,j)/L_th));
                    A(k-n,k) = mu_2; %west roosterpunt
                    A(k-1,k) = mu_2; %zuid roosterpunt
                    A(k+1,k) = mu_2; %noord roosterpunt
                elseif edge == "oost"
                    A(k,k) = -3*mu_2 + xi_2/C(i,j) * rho_0/rho_ster * (3*C(i,j) - C(i-1,j) - C(i,j+1)  - C(i+1,j)) - ((l^2/4)*(rho_in * L(i,j))/(1+L(i,j)/L_th));
                    A(k+n,k) = mu_2; %oost roosterpunt
                    A(k-1,k) = mu_2; %zuid roosterpunt
                    A(k+1,k) = mu_2; %noord roosterpunt
                elseif edge == "noord"
                    A(k,k) = -3*mu_2 + xi_2/C(i,j) * rho_0/rho_ster * (3*C(i,j) - C(i,j+1) - C(i,j-1) - C(i+1,j)) - ((l^2/4)*(rho_in * L(i,j))/(1+L(i,j)/L_th));
                    A(k+n,k) = mu_2; %oost roosterpunt
                    A(k-n,k) = mu_2; %west roosterpunt
                    A(k+1,k) = mu_2; %noord roosterpunt
                elseif edge == "zuid"
                    A(k,k) = -3*mu_2 + xi_2/C(i,j) * rho_0/rho_ster * (3*C(i,j) - C(i,j-1) - C(i,j+1) - C(i-1,j)) - ((l^2/4)*(rho_in * L(i,j))/(1+L(i,j)/L_th));
                    A(k+n,k) = mu_2; %oost roosterpunt
                    A(k-n,k) = mu_2; %west roosterpunt
                    A(k-1,k) = mu_2; %zuid roosterpunt
                elseif edge == "west_noord"
                    A(k,k) = -2*mu_2 + xi_2/C(i,j) * rho_0/rho_ster * (2*C(i,j) - C(i,j-1) - C(i+1,j)) - ((l^2/4)*(rho_in * L(i,j))/(1+L(i,j)/L_th));
                    A(k-n,k) = mu_2; %west roosterpunt
                    A(k+1,k) = mu_2; %noord roosterpunt
                elseif edge == "oost_noord"
                    A(k,k) = -2*mu_2 + xi_2/C(i,j) * rho_0/rho_ster * (2*C(i,j) - C(i+1,j) - C(i,j+1))- ((l^2/4)*(rho_in * L(i,j))/(1+L(i,j)/L_th));
                    A(k+n,k) = mu_2; %oost roosterpunt
                    A(k+1,k) = mu_2; %noord roosterpunt
                elseif edge == "west_zuid"
                    A(k,k) = -2*mu_2 + xi_2/C(i,j) * rho_0/rho_ster * (2*C(i,j) - C(i,j-1) - C(i-1,j)) - ((l^2/4)*(rho_in * L(i,j))/(1+L(i,j)/L_th));
                    A(k-n,k) = mu_2; %west roosterpunt
                    A(k-1,k) = mu_2; %zuid roosterpunt
                elseif edge == "oost_zuid"
                    A(k,k) = -2*mu_2 + xi_2/C(i,j) * rho_0/rho_ster * (2*C(i,j) - C(i-1,j) - C(i,j+1)) - ((l^2/4)*(rho_in * L(i,j))/(1+L(i,j)/L_th));
                    A(k+n,k) = mu_2; %oost roosterpunt
                    A(k-1,k) = mu_2; %zuid roosterpunt
                end
            end
        end
    end

    %% Vind de oplossing voor Au=f
    fs = A \ f;fs(isinf(fs))=0;%fs = reshape(fs,[n n]);
    output = reshape(f,[n n]);
end