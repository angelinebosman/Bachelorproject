function [output] = Macrofagen(n,r1,r2,m,C,L)
% Parameters
    mu_2 = 0.0005;
    xi_2 = 10;
    rho_0 = 5*10^(-7);
    rho_ster = 1;
    rho_m = 7;
    rho_in = 1*10^(-2);
    L_th = 1;
    


    l = 2*r1/n; %lengte tussen roosterpunten: diameter domein gedeeld door aantal roosterpunten
    A = sparse(n^2,n^2); %matrix die alleen de niet-nul waardes opslaat
    f = sparse(n^2,1);

    %discretisatie matrix
    for j= 1:n
        for i=1:n 
            k = i+(j-1)*n; %horizontale nummering
            [placement, edge] = indices(i,j,l,r1,r2); %bepaald waar in het domein ligt het roosterpunt
            if placement ~= "outside"
                f(k,1) = -(l^2/4 * rho_m * m(i,j));
            end
            if placement == "inside"
                A(k,k) = -3*mu_2 + (xi_2/C(i,j))*(rho_0/rho_ster)*(3*C(i,j) - C(i-1,j) - C(i,j-1) - C(i,j+1)) - (l^2/4)*(rho_in * L(i,j))/(1+L(i,j)/L_th);
                A(k+n,k) = mu_2; %oost roosterpunt
                A(k-n,k) = mu_2; %west roosterpunt
                A(k-1,k) = mu_2; %zuid roosterpunt
            elseif placement == "buitenrand"
                if edge == "west"
                    A(k,k) = -3*mu_2 + (xi_2/C(i,j))*(rho_0/rho_ster)*(2*C(i,j) - C(i-1,j) - C(i,j+1))- (l^2/4)*(rho_in * L(i,j))/(1+L(i,j)/L_th);
                    A(k+n,k) = mu_2; %oost roosterpunt
                    A(k-1,k) = mu_2; %zuid roosterpunt
                elseif edge == "oost"
                    A(k,k) = -3*mu_2 + (xi_2/C(i,j))*(rho_0/rho_ster)*(2*C(i,j) - C(i,j-1) - C(i-1,j))- (l^2/4)*(rho_in * L(i,j))/(1+L(i,j)/L_th);
                    A(k-n,k) = mu_2; %west roosterpunt
                    A(k-1,k) = mu_2; %zuid roosterpunt
                elseif edge == "noord"
                    A(k,k) = -3*mu_2 + (xi_2/C(i,j)) * rho_0/rho_ster * (3*C(i,j) - C(i-1,j) - C(i,j-1) - C(i,j+1)) -(l^2/4)*(rho_in * L(i,j))/(1+L(i,j)/L_th);
                    A(k-n,k) = mu_2; %west roosterpunt
                    A(k+n,k) = mu_2; %oost roosterpunt
                elseif edge == "zuid" 
                    A(k,k) = -3*mu_2 + (xi_2/C(i,j)) * rho_0/rho_ster * (2*C(i,j) - C(i,j+1) - C(i,j-1))- (l^2/4)*(rho_in * L(i,j))/(1+L(i,j)/L_th);
                    A(k-n,k) = mu_2; %west roosterpunt
                    A(k+n,k) = mu_2; %oost roosterpunt
                elseif edge == "oost_zuid"
                    A(k,k) = -2*mu_2 + (xi_2/C(i,j)) * rho_0/rho_ster * (C(i,j) - C(i,j-1)) - (l^2/4)*(rho_in * L(i,j))/(1+L(i,j)/L_th);
                    A(k-n,k) = mu_2; %west roosterpunt
                elseif edge == "west_zuid"
                    A(k,k) = -2*mu_2 + (xi_2/C(i,j)) * rho_0/rho_ster * (C(i,j) - C(i,j+1)) - (l^2/4)*(rho_in * L(i,j))/(1+L(i,j)/L_th);
                    A(k+n,k) = mu_2; %oost roosterpunt
                elseif edge == "west_noord"
                    A(k,k) =  -2*mu_2 + (xi_2/C(i,j)) * rho_0/rho_ster * (C(i,j) - C(i,j+1)) - (l^2/4)*(rho_in * L(i,j))/(1+L(i,j)/L_th);
                    A(k+n,k) = mu_2; %oost roosterpunt
                elseif edge == "oost_noord"
                    A(k,k) =  -2*mu_2 + (xi_2/C(i,j)) * rho_0/rho_ster * (2*C(i,j) - C(i-1,j) - C(i,j-1)) - (l^2/4)*(rho_in * L(i,j))/(1+L(i,j)/L_th);
                    A(k-n,k) = mu_2; %west roosterpunt
                    A(k-1,k) = mu_2; %zuid roosterpunt
                end
            elseif placement == "binnenrand"
                if edge == "west"
                    A(k,k) = -3*mu_2 + (xi_2/C(i,j)) * rho_0/rho_ster * (2*C(i,j) - C(i,j-1) - C(i-1,j)) - (l^2/4)*(rho_in * L(i,j))/(1+L(i,j)/L_th);
                    A(k-n,k) = mu_2; %west roosterpunt
                    A(k-1,k) = mu_2; %zuid roosterpunt
                elseif edge == "oost"
                    A(k,k) = -3*mu_2 + xi_2/C(i,j) * rho_0/rho_ster * (2*C(i,j) - C(i-1,j) - C(i,j+1)) - l^2/4*  (rho_in * L(i,j))/(1+L(i,j)/L_th);
                    A(k+n,k) = mu_2; %oost roosterpunt
                    A(k-1,k) = mu_2; %zuid roosterpunt
                elseif edge == "noord"
                    A(k,k) = -3*mu_2 + xi_2/C(i,j) * rho_0/rho_ster * (2*C(i,j) - C(i,j+1) - C(i,j-1)) - l^2/4*  (rho_in * L(i,j))/(1+L(i,j)/L_th);
                    A(k+n,k) = mu_2; %oost roosterpunt
                    A(k-n,k) = mu_2; %west roosterpunt
                elseif edge == "zuid"
                    A(k,k) = -3*mu_2 + xi_2/C(i,j) * rho_0/rho_ster * (3*C(i,j) - C(i-1,j) - C(i,j-1) - C(i,j+1)) - l^2/4*  (rho_in * L(i,j))/(1+L(i,j)/L_th);
                    A(k+n,k) = mu_2; %oost roosterpunt
                    A(k-n,k) = mu_2; %west roosterpunt
                elseif edge == "west_noord"
                    A(k,k) = -2*mu_2 + xi_2/C(i,j) * rho_0/rho_ster * (2*C(i,j) - C(i,j-1) - C(i-1,j)) - l^2/4*  (rho_in * L(i,j))/(1+L(i,j)/L_th);
                    A(k-n,k) = mu_2; %west roosterpunt
                    A(k-1,k) = mu_2; %zuid roosterpunt
                elseif edge == "oost_noord"
                    A(k,k) = -2*mu_2 + xi_2/C(i,j) * rho_0/rho_ster * (2*C(i,j) - C(i-1,j) - C(i,j+1)) - l^2/4*  (rho_in * L(i,j))/(1+L(i,j)/L_th);
                    A(k+n,k) = mu_2; %oost roosterpunt
                    A(k-1,k) = mu_2; %zuid roosterpunt
                elseif edge == "west_zuid"
                    A(k,k) = -2*mu_2 + xi_2/C(i,j) * rho_0/rho_ster * (2*C(i,j) - C(i,j-1) - C(i-1,j)) - l^2/4*  (rho_in * L(i,j))/(1+L(i,j)/L_th);
                    A(k-n,k) = mu_2; %west roosterpunt
                    A(k-1,k) = mu_2; %zuid roosterpunt
                elseif edge == "oost_zuid"
                    A(k,k) = -2*mu_2 + xi_2/C(i,j) * rho_0/rho_ster * (2*C(i,j) - C(i-1,j) - C(i,j+1)) - l^2/4*  (rho_in * L(i,j))/(1+L(i,j)/L_th);
                    A(k+n,k) = mu_2; %oost roosterpunt
                    A(k-1,k) = mu_2; %zuid roosterpunt
                end
            end
        end
    end


    fs = A \ f;fs(isinf(fs))=0;
    fs_reshape = reshape(fs,[n n]);
    output = fs_reshape
end