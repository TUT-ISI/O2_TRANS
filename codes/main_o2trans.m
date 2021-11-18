% ##############################################################################################
% #    Cooking recipe for oxygen transmittance computation                                     #
% ##############################################################################################
% This collection of scripts is developed to compute the oxygen
% transmittance spectrum in a defined optical path measured in [m]
% To compute the O2 transmittance, temperature(T) and pressure (p) measurements
% at a certain z. Always is required at minimum to have one measurement of
% (z,p,T) corresponding to the top of the optical path.
% 
% 
%   ---------------------------------- STEPS  ----------------------------------
% (1) Given three columns Z(m) p(atm) and T (K) I compute the air density
% (2) With the air density (rho_air) and the oxygen mixing ratio; compute the corresponding
% oxygen molecular volume density (O2rho)
% (3) Given (p,T) conditions, compute the oxygen transmittance for the assumed optical path
% (cos tetha corrected) 
% ----------------------------------------------------------------------------------
%
% This complement the Technical Note published in 
%
% Author is not responsible for any mistake or error in the computation.
% This is provided freely to the scientific comunity and for non-comercial use.
% It this scrips have been useful for your research work. Please, consider
% the citation. Thanks.
%
%
% Author: Neus Sabater
% Version v.0: 
% Version v1.0: Added Rayleigh transmittance, added CIA effect on O2-A
% Data: April/2020
% e-mail: neus.sabater@fmi.fi
% ----------------------------------------------------------------------------------
% ----------------------------------------------------------------------------------



% ---------------------
% INPUTS [*.txt file comma separated]
% --------------------
[num_samples,INPUT,GEO,path] = read_inputs;


% ---------------------
% Configure the settings
% ---------------------
flags = configure_settings();

f = waitbar(0,['Computing sample 1 of ', num2str(num_samples)]);

for i=1:num_samples   
    waitbar(i/num_samples,f,['Computing sample ',num2str(i),' of ',num2str(num_samples)]);
    % ---------------------
    % (1) Air density
    % ---------------------
    [rho_air,z_mes,p_mes,T_mes,zenith]= air_density(INPUT,GEO,i);
    T_ALL(i,:) = T_mes;
    
    
    % ---------------------
    % (2) O2 volume molecular density
    % ---------------------
    [O2rho]= O2_density(rho_air,z_mes,flags);
    rho_O2(i,:)= O2rho;
    
    % ---------------------
    % (3) O2 transmittance
    % ---------------------
    flags.iso = 16;
    [wl_O2B,wl_O2A,TRAN_O2B_O16,TRAN_O2A_O16,TRAN_O2B_RAY_O16,TRAN_O2A_RAY_O16] = O2_optical_depth(flags,O2rho,z_mes,p_mes,T_mes,zenith);
    flags.iso = 18;
    [wl_O2B_O16O18,wl_O2A_O16O18,TRAN_O2B_O16O18,TRAN_O2A_O16O18,TRAN_O2B_RAY_O16O18,TRAN_O2A_RAY_O16O18] = O2_optical_depth(flags,O2rho,z_mes,p_mes,T_mes,zenith);
    flags.iso = 17;
    [wl_O2B_O16O17,wl_O2A_O16O17,TRAN_O2B_O16O17,TRAN_O2A_O16O17,TRAN_O2B_RAY_O16O17,TRAN_O2A_RAY_O16O17] = O2_optical_depth(flags,O2rho,z_mes,p_mes,T_mes,zenith);

    
    % ---------------------
    % (4) Compute the CIA O2-Air transmittance
    % ---------------------    
    [t_cia_O2A,wl_cia_O2A,t_cia_O2B,wl_cia_O2B]= read_cia_file_hitran(flags.AUX,rho_air,O2rho,z_mes,T_mes); 
    TRAN_O2A_CIA   = interp1(wl_cia_O2A,t_cia_O2A',wl_O2A);   TRAN_O2A_CIA(isnan(TRAN_O2A_CIA))=1; 
    TRAN_O2B_CIA   = interp1(wl_cia_O2B,t_cia_O2B',wl_O2B);   TRAN_O2B_CIA(isnan(TRAN_O2B_CIA))=1; 
    
    % ---------------------
    % (5) Combine the total transmittance
    % ---------------------    
    TRAN_O2A_O16O18 = interp1(wl_O2A_O16O18,TRAN_O2A_O16O18,wl_O2A);   TRAN_O2A_O16O18(isnan(TRAN_O2A_O16O18))=1; 
    TRAN_O2A_O16O17 = interp1(wl_O2A_O16O17,TRAN_O2A_O16O17,wl_O2A);   TRAN_O2A_O16O17(isnan(TRAN_O2A_O16O17))=1; 
    T_O2A(i,:)      = TRAN_O2A_O16.*TRAN_O2A_O16O18.*TRAN_O2A_O16O17.*TRAN_O2A_CIA;
    
    TRAN_O2B_O16O18 = interp1(wl_O2B_O16O18,TRAN_O2B_O16O18,wl_O2B);   TRAN_O2B_O16O18(isnan(TRAN_O2B_O16O18))=1;
    TRAN_O2B_O16O17 = interp1(wl_O2B_O16O17,TRAN_O2B_O16O17,wl_O2B);   TRAN_O2B_O16O17(isnan(TRAN_O2B_O16O17))=1;
    T_O2B(i,:)      = TRAN_O2B_O16.*TRAN_O2B_O16O18.*TRAN_O2B_O16O17.*TRAN_O2B_CIA;
    
    TRAN_O2A_RAY_O16O18 = interp1(wl_O2A_O16O18,TRAN_O2A_RAY_O16O18,wl_O2A); TRAN_O2A_RAY_O16O18(isnan(TRAN_O2A_RAY_O16O18))=1;
    TRAN_O2A_RAY_O16O17 = interp1(wl_O2A_O16O17,TRAN_O2A_RAY_O16O17,wl_O2A); TRAN_O2A_RAY_O16O17(isnan(TRAN_O2A_RAY_O16O17))=1;
    T_O2A_RAY(i,:) = TRAN_O2A_RAY_O16.*TRAN_O2A_RAY_O16O18.*TRAN_O2A_RAY_O16O17;
    
    TRAN_O2B_RAY_O16O18 = interp1(wl_O2B_O16O18,TRAN_O2B_RAY_O16O18,wl_O2B); TRAN_O2B_RAY_O16O18(isnan(TRAN_O2B_RAY_O16O18))=1;
    TRAN_O2B_RAY_O16O17 = interp1(wl_O2B_O16O17,TRAN_O2B_RAY_O16O17,wl_O2B); TRAN_O2B_RAY_O16O17(isnan(TRAN_O2B_RAY_O16O17))=1;
    T_O2B_RAY(i,:) = TRAN_O2B_RAY_O16.*TRAN_O2B_RAY_O16O18.*TRAN_O2B_RAY_O16O17;
    
end
close(f)
f = waitbar(1,'Process finished :)');
close(f)
delete(f)
dlmwrite([path,'Transmittance_o2A.txt'],[wl_O2A',T_O2A'],'precision','%.6f');
dlmwrite([path,'Transmittance_o2B.txt'],[wl_O2B',T_O2B'],'precision','%.6f');
dlmwrite([path,'Transmittance_o2A_rayleigh.txt'],[wl_O2A',T_O2A_RAY'],'precision','%.16f');
dlmwrite([path,'Transmittance_o2B_rayleigh.txt'],[wl_O2B',T_O2B_RAY'],'precision','%.16f');

clear all
clc;


