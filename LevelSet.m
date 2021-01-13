function [output] = LevelSet(n,r1,dt,L,m,M,phi)
%% Assign paramater values
    mu_1 = 0.05;
    d1 = 0.1;
    sigma = 1;
    q = 7.5;
    m_th = 1;
    d2 = 0.005;
    mu_4 = 0.05;
    l = 2*r1/n;
    
 %% create fs as phi at timestep n+1
    fs = zeros(n,n);
    for i=2:(n-3)
        for j=2:(n-3)
            x = l*(i-1/2); %current position on the x-axis
            y = l*(j-1/2); %current position on the y-axis
            S = mu_1*(Lambda(sigma,L(i-1,j-1)) - d1*m(i-1,j-1)) + mu_4*(q*m(i-1,j-1)/(1+m(i-1,j-1) / m_th)*M(i-1,j-1) - d2*L(i-1,j-1));
            nabla_phi = 0;
            %% west en oost
            if x < 1/2
                nabla_phi = nabla_phi + phi(i,j) - phi(i-1,j); 
            elseif x > 3/2
                nabla_phi = nabla_phi + phi(i,j) - phi(i+1,j);
            end
             %% zuid en noord
                if y < 1/2
                    nabla_phi = nabla_phi + phi(i,j) - phi(i,j-1); 
                elseif y > 3/2
                    nabla_phi = nabla_phi + phi(i,j) - phi(i,j+1); 
                end
            v = -S;
            fs(i,j) = phi(i,j) + 3/(l^2)*v*dt*nabla_phi;
        end
    end

    output = fs;
end