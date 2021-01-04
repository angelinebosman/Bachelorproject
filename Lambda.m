function [output] = Lambda(sigma, LDL)
    %parameters
    gamma_0 = 20;
    sigma_0 = 2;
    L_th = 2;
    
    output = (gamma_0/(1+sigma/sigma_0))*(LDL/(1+LDL/L_th));
end