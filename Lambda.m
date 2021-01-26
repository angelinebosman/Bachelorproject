function [output] = Lambda(sigma, LDL)
    %% parameters
    gamma_0 = 2*10^4; %Snelheidsconstante waarmee de monocyten de intima betreden
    sigma_0 = 1; %nul shear stress snelheid van het endothelium
    L_th = 1; %Grenswaarde geoixideerd cholesterol
    
    %% Creëer output
    output = (gamma_0/(1+sigma/sigma_0))*(LDL/(1+LDL/L_th));
end