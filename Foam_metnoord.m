function [output] = Foam_metnoord(n,r1,r2,M,L,moving_bdry,xv, yv)
    if (~exist('xv', 'var')) %creates an optional input parameter
        opt = true;
    end
    if (~exist('yv', 'var'))
        opt = true;
    end
    
    %% Parameters
    mu_3 = 0.0005;
    rho_in = 10^(-2); 
    L_th = 1;

    %% create A and f
    l = 2*r1/n; %lengte tussen roosterpunten: diameter domein gedeeld door aantal roosterpunten
    A = sparse(n^2,n^2); %matrix die alleen de niet-nul waardes opslaat
    f = sparse(n^2,1);
    for j= 1:n
        for i=1:n
            k = i+(j-1)*n; %horizontale nummering
            if moving_bdry == 'y'
                [placement, edge] = indices_moving_bdry(i,j,l,r1, xv, yv);
            else
                [placement, edge] = indices(i,j,l,r1,r2); %bepaald waar in het domein het roosterpunt ligt
            end
            if placement ~= "outside"
                f(k,1) = -(l^2/4)*M(i,j)*(rho_in * L(i,j))/(1+L(i,j)/L_th);
            end
            if placement == "inside"
                A(k,k) = -4*mu_3;
                A(k+n,k) = mu_3; %oost roosterpunt
                A(k-n,k) = mu_3; %west roosterpunt
                A(k-1,k) = mu_3; %zuid roosterpunt
                A(k+1,k) = mu_3; %noord roosterpunt
            elseif placement == "buitenrand"
                if edge == "west"
                    A(k,k) = -3*mu_3;
                    A(k+n,k) = mu_3; %oost roosterpunt
                    A(k-1,k) = mu_3; %zuid roosterpunt
                    A(k+1,k) = mu_3; %noord roosterpunt
                elseif edge == "oost"
                    A(k,k) = -3*mu_3;
                    A(k-n,k) = mu_3; %west roosterpunt
                    A(k-1,k) = mu_3; %zuid roosterpunt
                    A(k+1,k) = mu_3; %noord roosterpunt
                elseif edge == "zuid"
                    A(k,k) = -3*mu_3;
                    A(k-n,k) = mu_3; %west roosterpunt
                    A(k+n,k) = mu_3; %oost roosterpunt
                    A(k+1,k) = mu_3; %noord roosterpunt
                elseif edge == "noord"
                    A(k,k) = -3*mu_3;
                    A(k-n,k) = mu_3; %west roosterpunt
                    A(k+n,k) = mu_3; %oost roosterpunt
                    A(k-1,k) = mu_3; %zuid roosterpunt
                elseif edge == "oost_zuid"
                    A(k,k) = -3*mu_3;
                    A(k-n,k) = mu_3; %west roosterpunt
                    A(k+1,k) = mu_3; %noord roosterpunt
                    A(k-1,k) = mu_3; %zuid roosterpunt
                elseif edge == "west_zuid"
                    A(k,k) = -3*mu_3;
                    A(k+n,k) = mu_3; %oost roosterpunt
                    A(k+1,k) = mu_3; %noord roosterpunt
                    A(k-1,k) = mu_3; %zuid roosterpunt
                elseif edge == "west_noord"
                    A(k,k) = -3*mu_3;
                    A(k+n,k) = mu_3; %oost roosterpunt
                    A(k-1,k) = mu_3; %zuid roosterpunt
                    A(k+1,k) = mu_3; %noord roosterpunt
                elseif edge == "oost_noord"
                    A(k,k) = -3*mu_3;
                    A(k-n,k) = mu_3; %west roosterpunt
                    A(k-1,k) = mu_3; %zuid roosterpunt
                    A(k+1,k) = mu_3; %noord roosterpunt
                end
            elseif placement == "binnenrand"
                if edge == "west"
                    A(k,k) = -3*mu_3; %binnenrand
                    A(k-n,k) = mu_3; %west roosterpunt
                    A(k-1,k) = mu_3; %zuid roosterpunt
                    A(k+1,k) = mu_3; %noord roosterpunt
                elseif edge == "oost"
                    A(k,k) = -3*mu_3; %binnenrand
                    A(k+n,k) = mu_3; %oost roosterpunt
                    A(k-1,k) = mu_3; %zuid roosterpunt
                    A(k+1,k) = mu_3; %noord roosterpunt
                elseif edge == "noord"
                    A(k,k) = -3*mu_3; %binnenrand
                    A(k+n,k) = mu_3; %oost roosterpunt
                    A(k-n,k) = mu_3; %west roosterpunt
                    A(k+1,k) = mu_3; %noord roosterpunt
                elseif edge == "zuid"
                    A(k,k) = -3*mu_3; %binnenrand
                    A(k+n,k) = mu_3; %oost roosterpunt
                    A(k-n,k) = mu_3; %west roosterpunt
                    A(k-1,k) = mu_3; %zuid roosterpunt
                elseif edge == "west_zuid"
                    A(k,k) = -3*mu_3; %binnenrand
                    A(k-n,k) = mu_3; %west roosterpunt
                    A(k-1,k) = mu_3; %zuid roosterpunt
                    A(k+1,k) = mu_3; %noord roosterpunt
                elseif edge == "oost_zuid"
                    A(k,k) = -3*mu_3; %binnenrand
                    A(k+n,k) = mu_3; %oost roosterpunt
                    A(k-1,k) = mu_3; %zuid roosterpunt
                    A(k+1,k) = mu_3; %noord roosterpunt
               elseif edge == "west_noord"
                    A(k,k) = -3*mu_3; %binnenrand
                    A(k-n,k) = mu_3; %west roosterpunt
                    A(k+1,k) = mu_3; %zuid roosterpunt
                    A(k+1,k) = mu_3; %noord roosterpunt
               elseif edge == "oost_noord"
                    A(k,k) = -3*mu_3; %binnenrand
                    A(k+n,k) = mu_3; %oost roosterpunt
                    A(k+1,k) = mu_3; %zuid roosterpunt
                    A(k+1,k) = mu_3; %noord roosterpunt
                end
            end
        end
    end

    %% find solution to Au=f
    fs = A \ f; fs(isinf(fs))=0;
    fs_reshape = reshape(fs, [n,n]);
    output = fs_reshape
end