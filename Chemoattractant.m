function [output] = Chemoattractant(n,r1,r2,m,M,F)
    % Parameters
    a1 = 0.5;
    a2 = 0.5;
    lamda = 30;
    mu_5 = 0.05;
    c0 = 1; %beginconcentratie


    l = 2*r1/n; %lengte tussen roosterpunten: diameter domein gedeeld door aantal roosterpunten
    A =sparse(n^2,n^2); %matrix die alleen de niet-nul waardes opslaat
    f = sparse(n^2,1);


    %discretisatie matrix
    for j= 1:n
        for i=1:n
            k = i+(j-1)*n; %horizontale nummering
            [placement, edge] = indices(i,j,l,r1,r2); %bepaald waar in het domein ligt het roosterpunt
            if placement ~= "outside"
                f(k,1) = -(l^2/4)*(lamda*F(i,j)) -c0;
            end
            if placement == "inside"
                A(k,k) = -3*mu_5 - (1/4*l^2)*(a1*m(i,j) + a2*M(i,j));
                A(k+n,k) = mu_5; %oost roosterpunt
                A(k-n,k) = mu_5; %west roosterpunt
                A(k-1,k) = mu_5; %zuid roosterpunt
            elseif placement == "buitenrand"
                if edge == "west"
                    A(k,k) = -3*mu_5- (1/4*l^2)*(a1*m(i,j) + a2*M(i,j));
                    A(k+n,k) = mu_5; %oost roosterpunt
                    A(k-1,k) = mu_5; %zuid roosterpunt
                elseif edge == "oost"
                    A(k,k) = -3*mu_5- (1/4*l^2)*(a1*m(i,j) + a2*M(i,j));
                    A(k-n,k) = mu_5; %west roosterpunt
                    A(k-1,k) = mu_5; %zuid roosterpunt
                elseif edge == "zuid"
                    A(k,k) = -3*mu_5- (1/4*l^2)*(a1*m(i,j) + a2*M(i,j));
                    A(k-n,k) = mu_5; %west roosterpunt
                    A(k+n,k) = mu_5; %oost roosterpunt
                elseif edge == "noord"
                    A(k,k) = -3*mu_5- (1/4*l^2)*(a1*m(i,j) + a2*M(i,j));
                    A(k-n,k) = mu_5; %west roosterpunt
                    A(k+n,k) = mu_5; %oost roosterpunt
                    A(k-1,k) = mu_5; %zuid roosterpunt
                elseif edge == "oost_zuid"
                    A(k,k) = -2*mu_5- (1/4*l^2)*(a1*m(i,j) + a2*M(i,j));
                    A(k-n,k) = mu_5; %west roosterpunt
                elseif edge == "west_zuid"
                    A(k,k) = -2*mu_5- (1/4*l^2)*(a1*m(i,j) + a2*M(i,j));
                    A(k+n,k) = mu_5; %oost roosterpunt
                elseif edge == "west_noord"
                    A(k,k) = -2*mu_5- (1/4*l^2)*(a1*m(i,j) + a2*M(i,j));
                    A(k+n,k) = mu_5; %oost roosterpunt
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
                elseif edge == "west_noord"
                    A(k,k) = -2*mu_5- (1/4*l^2)*(a1*m(i,j) + a2*M(i,j)); %binnenrand
                    A(k-n,k) = mu_5; %west roosterpunt
                    A(k-1,k) = mu_5; %zuid roosterpunt
                elseif edge == "oost"
                    A(k,k) = -3*mu_5- (1/4*l^2)*(a1*m(i,j) + a2*M(i,j)); %binnenrand
                    A(k+n,k) = mu_5; %oost roosterpunt
                    A(k-1,k) = mu_5; %zuid roosterpunt
                elseif edge == "oost_noord"
                    A(k,k) = -2*mu_5- (1/4*l^2)*(a1*m(i,j) + a2*M(i,j)); %binnenrand
                    A(k+n,k) = mu_5; %oost roosterpunt
                    A(k-1,k) = mu_5; %zuid roosterpunt
                elseif edge == "noord"
                    A(k,k) = -3*mu_5- (1/4*l^2)*(a1*m(i,j) + a2*M(i,j)); %binnenrand
                    A(k+n,k) = mu_5; %oost roosterpunt
                    A(k-n,k) = mu_5; %west roosterpunt
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


    fs = A\f;

    fs_reshape = reshape(fs,[n n]);
    output = fs_reshape
end