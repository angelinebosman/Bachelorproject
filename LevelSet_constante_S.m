function [output] = LevelSet_constante_S(n,r1,dt,S,phi)
%% initialisatie parameters
    l = 2*r1/n;
    
 %% creëer phi op tijdstap n+1
    %first_half = grad(:,1:n/2);
    %second_half = grad(:,n/2+1:n);
    [FX,FY] = gradient(phi);
    first = FY(1:n/4,1:n);
    second_part1 = FX(n/4+1:3*n/4,1:n/2);second_part2 = FX(n/4+1:3*n/4,n/2+1:n);second = [second_part1, -second_part2];
    third = -FY(3*n/4+1:n,1:n);
    grad = [first ; second ; third];
    output = phi - 3/(l^2).*dt.*S.*grad;
end