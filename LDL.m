function [output] = LDL(n,r1,r2,m,M,b)
    % Parameters
    rho_L = 2.5;
    mu_4 = 0.005;
    d2= 0.05;
    q = 7.5*10^5;
    m_th = 1;

    l = 2*r1/n; %lengte tussen roosterpunten: diameter domein gedeeld door aantal roosterpunten
    A = sparse(n^2,n^2); %matrix die alleen de niet-nul waardes opslaat
    f = zeros(n^2,1);


    %discretisatie matrix
    for j= 1:n
        for i=1:n
            k = i+(j-1)*n; %horizontale nummering
            [placement, edge] = indices(i,j,l,r1,r2); %bepaald waar in het domein ligt het roosterpunt
            if placement == "inside"
                A(k,k) = -3*mu_4 - (l^2/4)*rho_L*M(i,j);
                A(k+n,k) = mu_4; %oost roosterpunt
                A(k-n,k) = mu_4; %west roosterpunt
                A(k-1,k) = mu_4; %zuid roosterpunt
            elseif placement == "buitenrand"
                if edge == "west"
                    A(k,k) = -3*mu_4 - (l^2/4)*rho_L*M(i,j);
                    A(k+n,k) = mu_4; %oost roosterpunt
                    A(k-1,k) = mu_4; %zuid roosterpunt
                elseif edge == "oost"
                    A(k,k) = -3*mu_4- (l^2/4)*rho_L*M(i,j);
                    A(k-n,k) = mu_4; %west roosterpunt
                    A(k-1,k) = mu_4; %zuid roosterpunt
                elseif edge == "noord"
                    A(k,k) = -3*mu_4- (l^2/4)*rho_L*M(i,j);
                    A(k-n,k) = mu_4; %west roosterpunt
                    A(k+n,k) = mu_4; %oost roosterpunt
                elseif edge == "zuid"
                    A(k,k) = -3*mu_4- (l^2/4)*rho_L*M(i,j);
                    A(k-n,k) = mu_4; %west roosterpunt
                    A(k+n,k) = mu_4; %oost roosterpunt
                elseif edge == "oost_zuid"
                    A(k,k) = -2*mu_4- (l^2/4)*rho_L*M(i,j);
                    A(k-n,k) = mu_4; %west roosterpunt
                    A(k-1,k) = mu_4; %zuid roosterpunt
                elseif edge == "west_zuid"
                    A(k,k) = -2*mu_4 - (l^2/4)*rho_L*M(i,j);
                    A(k+n,k) = mu_4; %oost roosterpunt
                    A(k-1,k) = mu_4; %zuid roosterpunt
                elseif edge == "west_noord"
                    A(k,k) = -2*mu_4 - 1/4*l^2*rho_L*M(i,j);
                    A(k+n,k) = mu_4; %oost roosterpunt
                    A(k-1,k) = mu_4; %zuid roosterpunt
                elseif edge == "oost_noord"
                    A(k,k) = -3*mu_4 - (l^2/4)*rho_L*M(i,j);
                    A(k-n,k) = mu_4; %west roosterpunt
                    A(k-1,k) = mu_4; %zuid roosterpunt
                end
            elseif placement == "binnenrand"
                f(k,1) = -(sqrt(3)/3)*q*m(i,j)/(1+m(i,j) / m_th)*M(i,j);
                if edge == "west"
                    A(k,k) = -3*mu_4 - d2- (l^2/4)*rho_L*M(i,j); %binnenrand
                    A(k-n,k) = mu_4; %west roosterpunt
                    A(k-1,k) = mu_4; %zuid roosterpunt
                elseif edge == "oost"
                    A(k,k) = -3*mu_4 - d2 - (l^2/4)*rho_L*M(i,j); %binnenrand
                    A(k+n,k) = mu_4; %oost roosterpunt
                    A(k-1,k) = mu_4; %zuid roosterpunt
                elseif edge == "noord"
                    A(k,k) = -3*mu_4 - d2- (l^2/4)*rho_L*M(i,j); %binnenrand
                    A(k+n,k) = mu_4; %oost roosterpunt
                    A(k-n,k) = mu_4; %west roosterpunt
                elseif edge == "zuid"
                    A(k,k) = -3*mu_4 - d2- (l^2/4)*rho_L*M(i,j); %binnenrand
                    A(k+n,k) = mu_4; %oost roosterpunt
                    A(k-n,k) = mu_4; %west roosterpunt
                elseif edge == "oost_noord"
                    A(k,k) = -2*mu_4 - d2 - (l^2/4)*rho_L*M(i,j); %binnenrand
                    A(k+n,k) = mu_4; %oost roosterpunt
                    A(k-1,k) = mu_4; %zuid roosterpunt
                elseif edge == "oost_zuid"
                    A(k,k) = -2*mu_4 - d2 - (l^2/4)*rho_L*M(i,j); %binnenrand
                    A(k+n,k) = mu_4; %oost roosterpunt
                    A(k-1,k) = mu_4; %zuid roosterpunt
                elseif edge == "west_noord"
                    A(k,k) = -2*mu_4 - d2 - (l^2/4)*rho_L*M(i,j); %binnenrand
                    A(k-n,k) = mu_4; %west roosterpunt
                    A(k-1,k) = mu_4; %zuid roosterpunt
                elseif edge == "west_zuid"
                    A(k,k) = -2*mu_4 - d2- (l^2/4)*rho_L*M(i,j); %binnenrand
                    A(k-n,k) = mu_4; %west roosterpunt
                    A(k-1,k) = mu_4; %zuid roosterpunt
                end
            end
        end
    end
    
    b = reshape(b,[n*n,1]);
    fs = A\(f-b); fs(isinf(fs))=0;

    fs_reshape = reshape(fs,[n n]);
    %for i=1:n
    %    for j=1:n
    %        placement = indices(i,j,l,r1,r2);
    %        if placement == "outside"
    %            fs_reshape(i,j) = 0;
    %        end
    %    end
    %end
    output = fs_reshape
end