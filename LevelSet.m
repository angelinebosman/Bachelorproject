function [output] = LevelSet(n,r1,dt,L,m,M,phi,shear)
%% initialisatie parameters
    mu_1 = 0.05;
    d1 = 0.1;
    q = 7.5;
    m_th = 1;
    d2 = 0.005;
    mu_4 = 0.05;
    l = 2*r1/n;
    
 %% creëer phi op tijdstap n+1
    S = zeros(n,n);
    for i=2:(n-1)
        for j=2:(n-1)
            S(i,j) = mu_1*(Lambda(shear(i-1,j-1),L(i-1,j-1)) - d1*m(i-1,j-1)) + mu_4*(q*m(i-1,j-1)/(1+m(i-1,j-1) / m_th)*M(i-1,j-1) - d2*L(i-1,j-1));
        end
    end
    grad = gradient(phi);first_half = grad(:,1:n/2);second_half = grad(:,n/2+1:n);grad = [first_half, - second_half];
    output = phi - 4/(l^2).*dt.*S.*grad;
end