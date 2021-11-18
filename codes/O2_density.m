function  [O2rho]= O2_density(rho_air,z_mes,flags)
% O2_DENSITY This function computes the O2 density taking into account the
% air density and the O2 volume mixing ratio.

% INPUT: 
%     rho_air: Air density on each z step in [cm-3]
%     z_mes: Elevation above the surface level in [m]

% OUTPUT: 
%     O2rho: Oxygen volume density in [molecules · cm-3]
% Author: Neus Sabater
% Version v.0
% Data: April/2020
% e-mail: neus.sabater@fmi.fi
% ----------------------------------------------------------------------------------

A     = load([flags.AUX,filesep,'mol_const_I']);

Z_xi  = flipud(A(:,1)).*1000; % Z in m
xi_o2 = flipud(A(:,8)).*1e-6; % 1 ppmv = 1e-6 mol/mol

xi_o2_mes = interp1(Z_xi,xi_o2,z_mes);
O2rho = rho_air.*xi_o2_mes; % In [molecules /cm3]


end 
