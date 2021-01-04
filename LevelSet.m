function [output] = LevelSet(n,r1,r2,dt,L,m,M,b)
    %b is phi^n, oftewel de binnenrand in de vorige tijdsstap
    mu_1 = 0.05;
    d1 = 0.1;
    sigma = 1;
    q = 7.5;
    m_th = 1;
    d2 = 0.05;
    mu_4 = 0.005;
    
    l = 2*r1/n;
    A = sparse(n^2,n^2);
    f = (l^2/4)*b;
    
    for i=1:n
        for j=1:n
            k = i+(j-1)*n;
            [xv,yv] = get_xv_and_yv(b);
            [placement, edge] = indices_moving_bdry(i,j,l,r1,xv,yv);
            if placement == "binnenrand"
                S = mu_1*(Lambda(sigma,L(i,j)) - d1*m(i,j)) + mu_4*(q*m(i,j)/(1+m(i,j) / m_th)*M(i,j) - d2*L(i,j));
                if edge == "west"
                    A(k,k) = 2*dt*S;
                    A(k-n,k) = -dt*S; %west roosterpunt
                    A(k-1,k) = -dt*S; %zuid roosterpunt
                elseif edge == "west_noord"
                    A(k,k) = dt*S;
                    A(k-n,k) = -dt*S; %west roosterpunt
                elseif edge == "oost"
                    A(k,k) = 2*dt*S;
                    A(k+n,k) = -dt*S; %oost roosterpunt
                    A(k-1,k) = -dt*S; %zuid roosterpunt
                elseif edge == "oost_noord"
                    A(k,k) = 2*dt*S;
                    A(k+n,k) = -dt*S; %oost roosterpunt
                    A(k-1,k) = -dt*S; %zuid roosterpunt
                elseif edge == "noord"
                    A(k,k) = 2*dt*S;
                    A(k+n,k) = -dt*S; %oost roosterpunt
                    A(k-n,k) = -dt*S; %west roosterpunt
                elseif edge == "zuid"
                    A(k,k) = 3*dt*S;
                    A(k+n,k) = -dt*S; %oost roosterpunt
                    A(k-n,k) = -dt*S; %west roosterpunt
                    A(k-1,k) = -dt*S; %zuid roosterpunt
                elseif edge == "west_zuid"
                    A(k,k) = 2*dt*S;
                    A(k-n,k) = -dt*S; %west roosterpunt
                    A(k-1,k) = -dt*S; %zuid roosterpunt
                elseif edge == "oost_zuid"
                    A(k,k) = 2*dt*S;
                    A(k+n,k) = -dt*S; %oost roosterpunt
                    A(k-1,k) = -dt*S; %zuid roosterpunt
                end
            end
        end
    end
    f = reshape(f,[n*n,1]);
    fs = A*f;
    fs = reshape(fs,[n,n]);

    
    for i=1:n
        for j=1:n
            %placement = indices(i,j,l,r1,r2);
            if fs(i,j) < 0 %zit geen rand 
                fs(i,j) = 0;
            end
            if fs(i,j) > 0  %zit wel rand 
                fs(i,j) = 1;
            end
           %if placement == "binnenrand"
           %     if fs(i,j) == 0 %zit wel rand
           %        fs(i,j) = 1;
           %     end
           %end
            %if placement == "buitenrand"
            %    fs(i,j) = 1; %blijft altijd op dezelfde plek
            %end
        end
    end
    fs(isinf(fs))=0;
    output = fs;
end