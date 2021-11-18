function  [t_cia_O2A,wl_cia_O2A,t_cia_O2B,wl_cia_O2B] = read_cia_file_hitran(AUX,air_rho,O2_rho,path_length,T_mes)
% READ_CIA_FILE_HITRAN: This function reads the Colision Induced absorption coefficient file and
% computes the associated transmittance at a certain wavelengt range 

% INPUT: 
%
%     cia_file: Colision Induced absorption coefficient file from
%     https://hitran.org/cia/. See all the referenced cited in that
%     website. This function is generally programmed but assumes that we
%     are computing the cia effect of the O2-Air pair. Due to that reason
%     we need air and O2 densities as input. Please, refer to Eq 4 and 5
%     from the paper: Update of the HITRAN collision-induced absorption
%     section. Icarus. Karman et al. 2019.
%
%     air_rho: [molecule/cm^2]
%     O2_rho:  [molecule/cm^2]
%     path_length:  [cm]

%
% OUTPUT: 
%     wl: wavelength vector in [nm]
%     t_cia: Transmittance associated to the cia file par [nm]
% 
% Author: Neus Sabater
% Version v.0
% Data: Oct/2021
% e-mail: neus.sabater@fmi.fi
% ----------------------------------------------------------------------------------

% ---------------------
% For the O2-A region
% ---------------------

cia_T_files = [276:10:326]; % This info is hard coded. If Temperature ranges associated with the CIA files change these values should also change. 


for i=1:length(T_mes)-1 % Compute the CIA effects for each layer
    
    T_mes_mid        = T_mes(i) + ((T_mes(i+1)-T_mes(i))./2); % in K
    path_length_mid  = path_length(i)+ (path_length(i+1)-path_length(i))./2; % in m
    
    % Search the cia files required
    T_diff           = cia_T_files-T_mes(i);
    aux_T            = (T_diff<0) - (T_diff>=0);
    T_cia_index      = cia_T_files(diff(aux_T)~=0);
    
    if T_cia_index==326 % This is 50°C of air Tempea
        
        cia_file     = [AUX,filesep,'HITRAN',filesep,'O2-Aband_air_2021_T',num2str(T_cia_index),'K.cia'];
        M            = importdata(cia_file, ' ', 1);
        
        wv_num       = M.data(:,1); % [cm-1]
        k            = M.data(:,2); % [cm^5/molecule^2]   
    
    else
        cia_file     = [AUX,filesep,'HITRAN',filesep,'O2-Aband_air_2021_T',num2str(T_cia_index),'K.cia'];
        M            = importdata(cia_file, ' ', 1);     
        cia_file     = [AUX,filesep,'HITRAN',filesep,'O2-Aband_air_2021_T',num2str(T_cia_index+10),'K.cia'];
        M2           = importdata(cia_file, ' ', 1); 
        
        wv_num       = M.data(:,1); % [cm-1]
        k_1          = M.data(:,2); % [cm^5/molecule^2]   
        
        wv_num_2     = M2.data(:,1); % [cm-1]
        k_2          = M2.data(:,2); % [cm^5/molecule^2]   
        
        k_2          = interp1(wv_num_2,k_2,wv_num);            
        k            = interp1([T_cia_index, T_cia_index+10],[k_1';k_2'],T_mes_mid)';  
    end
    
  
    
    path_length_aux = path_length_mid.*100; % Convert from m to cm
        
    wl_cia_O2A = 1e7./wv_num;
    t_cia_O2A(i,:)  = exp(-(path_length_aux.*k.*air_rho(i).*O2_rho(i)));   % [-]    
end

clear k wv_num

% ---------------------
% For the O2-B region
% ---------------------

cia_T_files = [276,296]; % This info is hard coded. If Temperature ranges associated with the CIA files change these values should also change. 

for i=1:length(T_mes)-1 % Compute the CIA effects for each layer
    
        T_mes_mid        = T_mes(i) + ((T_mes(i+1)-T_mes(i))./2); % in K
        path_length_mid  = path_length(i)+ (path_length(i+1)-path_length(i))./2; % in m

        cia_file     = [AUX,filesep,'HITRAN',filesep,'O2-Bband_air_2021_T273K.cia'];
        M            = importdata(cia_file, ' ', 1);
        cia_file     = [AUX,filesep,'HITRAN',filesep,'O2-Bband_air_2021_T293K.cia'];
        M2           = importdata(cia_file, ' ', 1); 
        
        wv_num       = M.data(:,1); % [cm-1]
        k_1          = M.data(:,2); % [cm^5/molecule^2]   
        
        wv_num_2     = M2.data(:,1); % [cm-1]
        k_2          = M2.data(:,2); % [cm^5/molecule^2]   
        
        k_2          = interp1(wv_num_2,k_2,wv_num);      
        if T_mes(i)>=293            
            k            = interp1([273, 293],[k_1';k_2'],293)';
        else
            k            = interp1([273, 293],[k_1';k_2'],T_mes_mid)';
        end    
    
       path_length_aux = path_length_mid.*100; % Convert from m to cm
        
       wl_cia_O2B = 1e7./wv_num;
       t_cia_O2B(i,:)  = exp(-(path_length_aux.*k.*air_rho(i).*O2_rho(i)));   % [-]    
    
end








end