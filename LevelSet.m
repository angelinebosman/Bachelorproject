function [output] = LevelSet(n,r1,dt,L,m,M,phi)
%% initialisatie parameters
    %mu_1 = 0.05;
    d1 = 0.1;
    sigma = 1;
    q = 7.5;
    m_th = 1;
    d2 = 0.005;
    %mu_4 = 0.05;
    l = 2*r1/n;
    Size1 = 21;  %size macrofraag
    S_m = 16;  %size monocyte
    S_L = 2.4 / 1000; %size LDL particle
    
 %% creëer phi op tijdstap n+1
    S = zeros(n,n);
    for i=2:(n-1)
        for j=2:(n-1)
            %x = l*(i-1/2); %current position on the x-axis
            %y = l*(j-1/2); %current position on the y-axis
            S(i,j) = (S_L*(Lambda(sigma,L(i,j)) - d1*m(i,j)) + Size1*(q*m(i,j)/(1+m(i,j) / m_th)*M(i,j) - d2*L(i,j)))/(S_L + Size1);
            if S(i,j) < 0
                S(i,j) = abs(S(i,j));
            end
            %nabla_phi = 0;
            %if x < 1/2
            %    nabla_phi = nabla_phi + phi(i,j) - phi(i-1,j); 
            %elseif x > 1.25
            %    nabla_phi = nabla_phi + phi(i,j) - phi(i+1,j);
            %end
            %if y < 1/2
            %    nabla_phi = nabla_phi + phi(i,j) - phi(i,j-1); 
            %elseif y > 1.25
            %    nabla_phi = nabla_phi + phi(i,j) - phi(i,j+1); 
            %end
            %ps(i,j) = phi(i,j) - 3/(l^2)*S*dt*nabla_phi*10^90;
        end
    end
    grad = gradient(phi);first_half = grad(:,1:n/2);second_half = grad(:,n/2+1:n);grad = [first_half, - second_half];
    %[FX,FY] = gradient(phi);
    %first = FY(1:n/4,1:n);
    %second_part1 = FX(n/4+1:3*n/4,1:n/2);second_part2 = FX(n/4+1:3*n/4,n/2+1:n);second = [second_part1, -second_part2];
    %third = -FY(3*n/4+1:n,1:n);
    %grad = [first ; second ; third];
    output = phi - 3/(l^2).*dt.*S.*grad;
end