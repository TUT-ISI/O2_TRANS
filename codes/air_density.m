function  [rho_air,z_mes,p_mes,T_mes,zenith]= air_density(AUX,GEO,i)
% AIR_DENSITY This function computes the air density following the ideal
% gas law and according to the pressure(p) and temperature (T) values
% contained in the INPUT file

% INPUT: 
%     AUX amtrix must contain at least three columns being 
%     ----- sample 1 ------
%       1st column Z [m]
%       2nd column p [mbar]
%       3rd column T [K]
%     ----- sample 2 ------
%       4th column Z [m]
%       5th column p [mbar]
%       6th column T [K]
%     ----- sample i ------
%     ((i-1)*3)+1 column Z [m]
%     ((i-1)*3)+2 column p [mbar]
%     ((i-1)*3)+3 column T [K]
%     ------------------------
%     i the set of triplets evaluated; i.e. data samples

% OUTPUT: 
%     rho_air: Air density [cm^-3]
%     z_mes: elevation above the surface in [m]
%     p_mes: pressure at each z in [hPa or mbar]
%     T_mes: Temperature at each z in [K]
% Author: Neus Sabater
% Version v.0
% Data: April/2020
% e-mail: neus.sabater@fmi.fi
% ----------------------------------------------------------------------------------


z_mes = AUX(:,((i-1)*3)+1); % In [m]
p_mes = AUX(:,((i-1)*3)+2); % In [hPa]
T_mes = AUX(:,((i-1)*3)+3); % In [K]

g     = 9.8; % Correction with the latitude is needed?
R     = 8.3144e4; % [cm3 hPa/K mol]
Na    = 6.022e23; %[-]
K     = 1.3807e-19; % k = R./Na 

rho_air = p_mes./(K.*T_mes);
zenith  = GEO(i); %The Zenith angle in degree associated with that profile

end
