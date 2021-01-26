function [output] = LDL(n,r1,r2,m,M,b,homogeen,moving_bdry,xv, yv)
    if (~exist('xv', 'var')) %creates an optional input parameter
        opt = true;
    end
    if (~exist('yv', 'var'))
        opt = true;
    end
    %% Parameters
    if homogeen == 'y'
        rho_L = 2.5; %Verteringsconstante van geoixideerde LDL-deeltjes door macrofagen
    else
        rho_L = 250;%Verteringsconstante van geoixideerde LDL-deeltjes door macrofagen
    end
    mu_4 = 0.005; %chemotactische gevoeligheid
    d2= 0.125; %diffusieconstante van bloedwand naar lumen
    q = 7.5; %opnamesnelheid LDL-deeltjes door het endothelium
    m_th = 1;

    %% create A and f
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
                A(k,k) = -4*mu_4 - (l^2/4)*rho_L*M(i,j);
                A(k+n,k) = mu_4; %oost roosterpunt
                A(k-n,k) = mu_4; %west roosterpunt
                A(k-1,k) = mu_4; %zuid roosterpunt
                A(k+1,k) = mu_4; %noord roosterpunt
            elseif placement == "buitenrand"
                if edge == "west"
                    A(k,k) = -3*mu_4 - (l^2/4)*rho_L*M(i,j);
                    A(k+n,k) = mu_4; %oost roosterpunt
                    A(k-1,k) = mu_4; %zuid roosterpunt
                    A(k+1,k) = mu_4; %noord roosterpunt
                elseif edge == "oost"
                    A(k,k) = -3*mu_4- (l^2/4)*rho_L*M(i,j);
                    A(k-n,k) = mu_4; %west roosterpunt
                    A(k-1,k) = mu_4; %zuid roosterpunt
                    A(k+1,k) = mu_4; %noord roosterpunt
                elseif edge == "noord"
                    A(k,k) = -3*mu_4- (l^2/4)*rho_L*M(i,j);
                    A(k-n,k) = mu_4; %west roosterpunt
                    A(k+n,k) = mu_4; %oost roosterpunt
                    A(k-1,k) = mu_4; %zuid roosterpunt
                elseif edge == "zuid"
                    A(k,k) = -3*mu_4- (l^2/4)*rho_L*M(i,j);
                    A(k-n,k) = mu_4; %west roosterpunt
                    A(k+n,k) = mu_4; %oost roosterpunt
                    A(k+1,k) = mu_4; %noord roosterpunt
                elseif edge == "oost_zuid"
                    A(k,k) = -2*mu_4- (l^2/4)*rho_L*M(i,j);
                    A(k-n,k) = mu_4; %west roosterpunt
                    A(k+1,k) = mu_4; %noord roosterpunt
                elseif edge == "west_zuid"
                    A(k,k) = -2*mu_4 - (l^2/4)*rho_L*M(i,j);
                    A(k+n,k) = mu_4; %oost roosterpunt
                    A(k+1,k) = mu_4; %noord roosterpunt
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
                    A(k+1,k) = mu_4; %noord roosterpunt
                elseif edge == "oost"
                    A(k,k) = -3*mu_4 - d2 - (l^2/4)*rho_L*M(i,j); %binnenrand
                    A(k+n,k) = mu_4; %oost roosterpunt
                    A(k-1,k) = mu_4; %zuid roosterpunt
                    A(k+1,k) = mu_4; %noord roosterpunt
                elseif edge == "noord"
                    A(k,k) = -3*mu_4 - d2- (l^2/4)*rho_L*M(i,j); %binnenrand
                    A(k+n,k) = mu_4; %oost roosterpunt
                    A(k-n,k) = mu_4; %west roosterpunt
                    A(k+1,k) = mu_4; %noord roosterpunt
                elseif edge == "zuid"
                    A(k,k) = -3*mu_4 - d2- (l^2/4)*rho_L*M(i,j); %binnenrand
                    A(k+n,k) = mu_4; %oost roosterpunt
                    A(k-n,k) = mu_4; %west roosterpunt
                    A(k-1,k) = mu_4; %zuid roosterpunt
                elseif edge == "oost_noord"
                    A(k,k) = -2*mu_4 - d2 - (l^2/4)*rho_L*M(i,j); %binnenrand
                    A(k+n,k) = mu_4; %oost roosterpunt
                    A(k+1,k) = mu_4; %noord roosterpunt
                elseif edge == "oost_zuid"
                    A(k,k) = -2*mu_4 - d2 - (l^2/4)*rho_L*M(i,j); %binnenrand
                    A(k+n,k) = mu_4; %oost roosterpunt
                    A(k-1,k) = mu_4; %zuid roosterpunt
                elseif edge == "west_noord"
                    A(k,k) = -2*mu_4 - d2 - (l^2/4)*rho_L*M(i,j); %binnenrand
                    A(k-n,k) = mu_4; %west roosterpunt
                    A(k+1,k) = mu_4; %noord roosterpunt
                elseif edge == "west_zuid"
                    A(k,k) = -2*mu_4 - d2- (l^2/4)*rho_L*M(i,j); %binnenrand
                    A(k-n,k) = mu_4; %west roosterpunt
                    A(k-1,k) = mu_4; %zuid roosterpunt
                end
            end
        end
    end
    
    %% find solution to Au=f
    b = reshape(b,[n*n,1]);
    fs = A\(f-b);fs(isinf(fs))=0;
    fs_reshape = reshape(fs,[n, n]);
    output = abs(fs_reshape)
end