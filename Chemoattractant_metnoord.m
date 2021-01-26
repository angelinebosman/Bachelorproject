function [output] = Chemoattractant(n,r1,r2,m,M,F,moving_bdry,xv, yv)
    if (~exist('xv', 'var')) %creates an optional input parameter
        opt = true;
    end
    if (~exist('yv', 'var'))
        opt = true;
    end
    
    %% Parameters
    a1 = 0.5;
    a2 = 0.5;
    lamda = 30;
    mu_5 = 0.05;
    c0 = 10^(-10); %beginconcentratie
    
    %% create A and f
    l = 2*r1/n; %lengte tussen roosterpunten: diameter domein gedeeld door aantal roosterpunten
    A =sparse(n^2,n^2); %matrix die alleen de niet-nul waardes opslaat
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
                f(k,1) = -(l^2/4)*(lamda*F(i,j)) - c0;
            end
            if placement == "inside"
                A(k,k) = -4*mu_5 - (1/4*l^2)*(a1*m(i,j) + a2*M(i,j));
                A(k+n,k) = mu_5; %oost roosterpunt
                A(k-n,k) = mu_5; %west roosterpunt
                A(k-1,k) = mu_5; %zuid roosterpunt
                A(k+1,k) = mu_5; %noord roosterpunt
            elseif placement == "buitenrand"
                if edge == "west"
                    A(k,k) = -3*mu_5 - (1/4*l^2)*(a1*m(i,j) + a2*M(i,j));
                    A(k+n,k) = mu_5; %oost roosterpunt
                    A(k-1,k) = mu_5; %zuid roosterpunt
                    A(k+1,k) = mu_5; %noord roosterpunt
                elseif edge == "oost"
                    A(k,k) = -3*mu_5- (1/4*l^2)*(a1*m(i,j) + a2*M(i,j));
                    A(k-n,k) = mu_5; %west roosterpunt
                    A(k-1,k) = mu_5; %zuid roosterpunt
                    A(k+1,k) = mu_5; %noord roosterpunt
                elseif edge == "zuid"
                    A(k,k) = -3*mu_5- (1/4*l^2)*(a1*m(i,j) + a2*M(i,j));
                    A(k-n,k) = mu_5; %west roosterpunt
                    A(k+n,k) = mu_5; %oost roosterpunt
                    A(k+1,k) = mu_5; %noord roosterpunt
                elseif edge == "noord"
                    A(k,k) = -3*mu_5- (1/4*l^2)*(a1*m(i,j) + a2*M(i,j));
                    A(k-n,k) = mu_5; %west roosterpunt
                    A(k+n,k) = mu_5; %oost roosterpunt
                    A(k-1,k) = mu_5; %zuid roosterpunt
                elseif edge == "oost_zuid"
                    A(k,k) = -2*mu_5- (1/4*l^2)*(a1*m(i,j) + a2*M(i,j));
                    A(k-n,k) = mu_5; %west roosterpunt
                    A(k+1,k) = mu_5; %noord roosterpunt
                elseif edge == "west_zuid"
                    A(k,k) = -2*mu_5- (1/4*l^2)*(a1*m(i,j) + a2*M(i,j));
                    A(k+n,k) = mu_5; %oost roosterpunt
                    A(k+1,k) = mu_5; %noord roosterpunt
                elseif edge == "west_noord"
                    A(k,k) = -2*mu_5- (1/4*l^2)*(a1*m(i,j) + a2*M(i,j));
                    A(k+n,k) = mu_5; %oost roosterpunt
                    A(k-1,k) = mu_5; %zuid roosterpunt
                elseif edge == "oost_noord"
                    A(k,k) = -2*mu_5- (1/4*l^2)*(a1*m(i,j) + a2*M(i,j));
                    A(k-n,k) = mu_5; %west roosterpunt
                    A(k-1,k) = mu_5; %zuid roosterpunt
                end
            elseif placement == "binnenrand"
                if edge == "west"
                    A(k,k) = -3*mu_5- (1/4*l^2)*(a1*m(i,j) + a2*M(i,j)) ; %binnenrand
                    A(k-n,k) = mu_5; %west roosterpunt
                    A(k-1,k) = mu_5; %zuid roosterpunt
                    A(k+1,k) = mu_5; %noord roosterpunt
                elseif edge == "west_noord"
                    A(k,k) = -2*mu_5- (1/4*l^2)*(a1*m(i,j) + a2*M(i,j)); %binnenrand
                    A(k-n,k) = mu_5; %west roosterpunt
                    A(k+1,k) = mu_5; %noord roosterpunt
                elseif edge == "oost"
                    A(k,k) = -3*mu_5- (1/4*l^2)*(a1*m(i,j) + a2*M(i,j)); %binnenrand
                    A(k+n,k) = mu_5; %oost roosterpunt
                    A(k-1,k) = mu_5; %zuid roosterpunt
                    A(k+1,k) = mu_5; %noord roosterpunt
                elseif edge == "oost_noord"
                    A(k,k) = -2*mu_5- (1/4*l^2)*(a1*m(i,j) + a2*M(i,j)); %binnenrand
                    A(k+n,k) = mu_5; %oost roosterpunt
                    A(k+1,k) = mu_5; %noord roosterpunt
                elseif edge == "noord"
                    A(k,k) = -3*mu_5- (1/4*l^2)*(a1*m(i,j) + a2*M(i,j)); %binnenrand
                    A(k+n,k) = mu_5; %oost roosterpunt
                    A(k-n,k) = mu_5; %west roosterpunt
                    A(k+1,k) = mu_5; %noord roosterpunt
                elseif edge == "zuid"
                    A(k,k) = -3*mu_5- (1/4*l^2)*(a1*m(i,j) + a2*M(i,j)); %binnenrand
                    A(k+n,k) = mu_5; %oost roosterpunt
                    A(k-n,k) = mu_5; %west roosterpunt
                    A(k-1,k) = mu_5; %zuid roosterpunt
                elseif edge == "west_zuid"
                    A(k,k) = -2*mu_5- (1/4*l^2)*(a1*m(i,j) + a2*M(i,j)); %binnenrand
                    A(k-n,k) = mu_5; %west roosterpunt
                    A(k-1,k) = mu_5; %zuid roosterpunt
                elseif edge == "oost_zuid"
                    A(k,k) = -2*mu_5- (1/4*l^2)*(a1*m(i,j) + a2*M(i,j)); %binnenrand
                    A(k+n,k) = mu_5; %oost roosterpunt
                    A(k-1,k) = mu_5; %zuid roosterpunt
                end
            end
        end
    end

    %% find solution to Au=f
    fs = A\f;
    fs_reshape = reshape(fs,[n n]);
    fs_reshape(isinf(fs_reshape)) = 0;
    output = fs_reshape
end