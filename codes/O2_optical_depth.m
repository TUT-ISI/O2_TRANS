function [wl_O2B,wl_O2A,TRAN_O2B,TRAN_O2A,Total_tran_B_ray,Total_tran_A_ray] = O2_optical_depth(flags,O2rho,z_mes,p_mes,T_mes,zenith)
% O2_DENSITY This function computes the O2 density taking into account the
% air density and the O2 volume mixing ratio.

% INPUT: 
%     flags: flags containing the settings
%     O2rho: Oxygen volume density in [molecules · cm-3]
%     z_mes: Elevation above the surface level in [m]
%     p_mes: pressure at each z in [hPa or mbar]
%     T_mes: Temperature at each z in [K]
%
% OUTPUT: 
%     wl_O2B: wavelength vector in [nm]
%     wl_O2A: wavelength vector in [nm]
%     TRAN_O2B: oxygen transmittance in the O2-B region (corresponding to wl_O2B)
%     TRAN_O2A: oxygen transmittance in the O2-A region (corresponding to wl_O2A)
% 
% Author: Neus Sabater
% Version v.0
% Data: April/2020
% e-mail: neus.sabater@fmi.fi
% ----------------------------------------------------------------------------------


% ---------------------
% Internal settings
% ---------------------
figdir    = flags.dir;
flag_plot = flags.plotting; % [1]= YES, [other]=NO
flag_save = flags.saving;   % [1]= YES, [other]=NO
% ---------------------

% -----------------------------------------------------------
% (1) Read the *.par file from HITRAN
% -----------------------------------------------------------

if flags.iso == 16
    par_file  = [flags.AUX,filesep,'HITRAN',filesep,'616fccbd_O16O16.par']; % isotopologue 16O16O 
elseif flags.iso == 17
    par_file  = [flags.AUX,filesep,'HITRAN',filesep,'616fcc95_O16O17.par']; % isotopologue 16O17O 
elseif flags.iso == 18
    par_file  = [flags.AUX,filesep,'HITRAN',filesep,'616fcc34_O16O18.par']; % isotopologue 16O18O 
end
[lines]   = importhitran(par_file);

% -----------------------------------------------------------
% Define the constants needed in HITRAN (not IS)
% -----------------------------------------------------------
cons.h           = 6.62606957e-27; % [erg·s] Planck constant 
cons.c           = 2.99792458e10; % [cm/s] Speed of Light 
cons.k           = 1.3806488e-16; %[erg/K] Boltzmann constant 
cons.c2          = 1.4387770; % [cm · K] C2 = hc/k second radiation constant
cons.Na          = 6.02214129e23; %[mol-1] Avogadro number
cons.p_ref       = 1; % [atm]
cons.T_ref       = 296; % [K]
cons.g           = 9.8;% [m/s2] % This could depend on the lat
cons.co2         = 0.00036; % 360ppm 
cons.ma          = 15.0556.*cons.co2+28.9595; %g/mol Expression from On Rayleigh optical depth calculations. Barry Bodhaine.



% These depend on the isotopologe
if flags.iso == 16
    cons.abundance   = 9.95262E-01;
    cons.O2_mol_frac = 0.20946.*cons.abundance; % [mol/mol_air]
    cons.M           = 31.989830; % Molar mass of the isotopologue (O2-16) in [g]
    cons.Q           = 2.1573e2; % Q(T=296 K)

elseif flags.iso == 17
    cons.abundance   = 7.422350e-4;
    cons.O2_mol_frac = 0.20946.*cons.abundance; % [mol/mol_air]
    cons.M           = 32.994045; % Molar mass of the isotopologue (O2-16) in [g]
    cons.Q           = 2658.12; % Q(T=296 K) 2.6581e3
   
elseif flags.iso == 18
    cons.abundance   = 0.003991;
    cons.O2_mol_frac = 0.20946.*cons.abundance; % [mol/mol_air]
    cons.M           = 33.994076; % Molar mass of the isotopologue (O2-16) in [g]
    cons.Q           = 455.23; % Q(T=296 K) 4.5523e2

end
% The sum of all the abundaces gives 1


% -----------------------------------------------------------
% load partition sums
% -----------------------------------------------------------
%partt_file = [flags.AUX,filesep,'HITRAN',filesep,'q36.txt']; % O2 isotope 16

if flags.iso == 16
    partt_file = [pwd,filesep,'AUX',filesep,'HITRAN',filesep,'q36_O16.txt'];    % O2 isotopologue O16O16
elseif flags.iso == 17
    partt_file = [pwd,filesep,'AUX',filesep,'HITRAN',filesep,'q38_O16O17.txt']; % O2 isotopologue O16O17
elseif flags.iso == 18
    partt_file = [pwd,filesep,'AUX',filesep,'HITRAN',filesep,'q37_O16O18.txt']; % O2 isotopologue O16O18
end


aux_Q_T    = importdata (partt_file);
TT = aux_Q_T(:,1);
QQ = aux_Q_T(:,2);
% -----------------------------------------------------------
% Convert the wavenumbers (cm-1) into nm
% -----------------------------------------------------------
nm =1e7./lines.transitionWavenumber;

% -----------------------------------------------------------
% Find the limits for the O2-B and the O2-A bands
% -----------------------------------------------------------
[~,O2B_ini] = min(abs(nm-684)); 
[~,O2B_end] = min(abs(nm-698)); 
range_B = O2B_ini-O2B_end+1;

[~,O2A_ini] = min(abs(nm-758)); 
[~,O2A_end] = min(abs(nm-770)); 
range_A = O2A_ini-O2A_end+1;



if flag_plot
    % plotting the Intensity Line at T=296 K
    figure
    subplot(1,2,1)
    for i=O2B_end:O2B_ini
        plot([nm(i) nm(i)],[lines.lineIntensity(i),0],'r')
        hold on
    end
    xlabel('\lambda [nm]')
    ylabel('S_{ij}(T_{ref}) [cm^{-1}/(molecule cm^2)]')
    title('O_2-B region')
    set(gca,'Fontsize',18)
    
    subplot(1,2,2)
    for i=O2A_end:O2A_ini
        plot([nm(i) nm(i)],[lines.lineIntensity(i),0],'r')
        hold on
    end
    xlabel('\lambda [nm]')
    ylabel('S_{ij}(T_{ref}) [cm^{-1}/(molecule cm^2)]')
    title('O_2-A region')
    set(gca,'Fontsize',18)
    
    set(gcf,'PaperPositionMode','auto')
    set(gcf, 'Position', [300 500 1500 400])
        
        if flag_save            
            saveas(gcf, [figdir,filesep,'S_ij_Tref_O16O',num2str(flags.iso),'.fig'],'fig');
            saveas(gcf, [figdir,filesep,'S_ij_Tref_O16O',num2str(flags.iso),'.eps'],'epsc');
            saveas(gcf, [figdir,filesep,'S_ij_Tref_O16O',num2str(flags.iso),'.png'],'png');
        end
end


% -----------------------------------------------------------
% Plotting the approximation of the Lorentizian functions at a
% high spectral resolution for the T=296 K and p= 1atm
% Taking into account the line intensity Sij (integral) and the 
% air broadenned halfwidth T=296 K
% ----------------------------------------------------------
% All is computed in wavenumbers [cm-1] and later plotted in wavelengths [nm]
% Here we generate the High Resolution (HR) matrix of nm_hr_O2B, and nm_hr_O2A
% where dim 1 counts for the number of absorption lines and dim 2 counts
% for the points representing the absorption intensity function.
% -----Using the default values in the HAPI, 10 cm-1 wings-------
cont = 0;
for i = O2B_end:O2B_ini
    cont = cont+1;
    v_hr_O2B(cont,:) = linspace(lines.transitionWavenumber(i)-10,lines.transitionWavenumber(i)+10,3000);
end
nm_hr_O2B = 1e7./v_hr_O2B;

cont = 0;
for i = O2A_end:O2A_ini
    cont = cont+1;
    v_hr_O2A(cont,:) = linspace(lines.transitionWavenumber(i)-10,lines.transitionWavenumber(i)+10,3000);
end
nm_hr_O2A = 1e7./v_hr_O2A;
% -----------------------



cont = 0;
for i = O2B_end:O2B_ini
    cont = cont+1;
    L_B(cont,:)      = lines.lineIntensity(i).*(lines.airBroadenedWidth(i)./(pi*((v_hr_O2B(cont,:)-lines.transitionWavenumber(i)).^2 + lines.airBroadenedWidth(i).^2)));
    L_B_self(cont,:) = lines.lineIntensity(i).*(lines.selfBroadenedWidth(i)./(pi*((v_hr_O2B(cont,:)-lines.transitionWavenumber(i)).^2 + lines.selfBroadenedWidth(i).^2)));
end


%  ------
% | O2-B |
%  ------




if flag_plot
    figure      % Air boradened
    for cont=1:range_B
        plot(nm_hr_O2B(cont,:),L_B(cont,:),'color','b')
        hold on
    end
    for i=O2B_end:O2B_ini
        plot([nm(i) nm(i)],[lines.lineIntensity(i)./(pi.*lines.airBroadenedWidth(i)),0],'r') % Maximun for a Lorentzian function (1/pi*gamma_HI)
        hold on
    end
    
    for cont=1:range_B
        plot(nm_hr_O2B(cont,:),L_B_self(cont,:),'color', 'k')
        hold on
    end
    grid on
    xlabel('\lambda [nm]')
    ylabel('Air boradened [cm^{-1}/atm]')
    title('O_2-B air Broadened')
    set(gca,'Fontsize',18)
    set(gcf,'PaperPositionMode','auto')
    set(gcf, 'Position', [300 500 1500 400])
    if flag_save
        saveas(gcf, [figdir,filesep,'Air_broad_B_O16O',num2str(flags.iso),'.fig'],'fig');
        saveas(gcf, [figdir,filesep,'Air_broad_B_O16O',num2str(flags.iso),'.eps'],'epsc');
        saveas(gcf, [figdir,filesep,'Air_broad_B_O16O',num2str(flags.iso),'.png'],'png');
    end
end




%  ------
% | O2-A |
%  ------

cont =0;
for i = O2A_end:O2A_ini
     cont = cont+1;
    L_A(cont,:)      = lines.lineIntensity(i).*(lines.airBroadenedWidth(i)./(pi*((v_hr_O2A(cont,:)-lines.transitionWavenumber(i)).^2 + lines.airBroadenedWidth(i).^2)));
    L_A_self(cont,:) = lines.lineIntensity(i).*(lines.selfBroadenedWidth(i)./(pi*((v_hr_O2A(cont,:)-lines.transitionWavenumber(i)).^2 + lines.selfBroadenedWidth(i).^2)));
end


if flag_plot
    figure    
    for cont=1:range_A
        plot(nm_hr_O2A(cont,:),L_A(cont,:),'color', 'b')
        hold on
    end
    for i=O2A_end:O2A_ini
        plot([nm(i) nm(i)],[lines.lineIntensity(i)./(pi.*lines.airBroadenedWidth(i)),0],'r') % Maximun for a Lorentzian function (1/pi*gamma_HI)
        hold on
    end 
    
    set(gca,'Fontsize',18)
    for cont=1:range_A
        plot(nm_hr_O2A(cont,:),L_A_self(cont,:),'color','k')
        hold on
    end
    set(gca,'Fontsize',18)
    grid on
    xlabel('\lambda [nm]')
    ylabel('Air boradened [cm^{-1}/atm]')
    title('O_2-A air Broadened')
    set(gcf,'PaperPositionMode','auto')
    set(gcf, 'Position', [300 500 1500 400])
    
    if flag_save
        saveas(gcf, [figdir,filesep,'Air_broad_B_O16O',num2str(flags.iso),'.fig'],'fig');
        saveas(gcf, [figdir,filesep,'Air_broad_B_O16O',num2str(flags.iso),'.eps'],'epsc');
        saveas(gcf, [figdir,filesep,'Air_broad_B_O16O',num2str(flags.iso),'.png'],'png');
    end
end


% -----------------------------------------------------------
% TEMPERATURE and PRESSURE dependence of the line width
% Computations need to be done in cm-1. Later we can plot in [nm]
% -----------------------------------------------------------

% Compute the average Temperature and pressure on each interval measured 
if length(T_mes)==1    % If there is only one measurement point (not a profile)
    T      = T_mes; %[K]
    p      = p_mes; % [mbar]
    p      = p./1013.25; % from mbar to atm
else  % If there is a profile   
    T      = T_mes(1:end-1)+(diff(T_mes)./2); %[K]
    p      = p_mes(1:end-1)+(diff(p_mes)./2); % [mbar]
    p      = p./1013.25; % from mbar to atm
end


Q_T    = interp1(TT,QQ,T);
% p      = cons.p_ref.*T./cons.T_ref; % Isochoric process P1*v1 = NR*T1;P2*v2 = NR*T2. [atm]
p_self = cons.O2_mol_frac.*p;

cont = 0;
for i=O2B_end:O2B_ini
    cont = cont+1;
    alpha_D_O2B(cont,:) = (lines.transitionWavenumber(i)/(cons.c)).*sqrt((2.*cons.Na.*cons.k.*T.*log(2))./(cons.M));
end

cont = 0;
for i=O2A_end:O2A_ini
    cont = cont+1;
    alpha_D_O2A(cont,:) = (lines.transitionWavenumber(i)/(cons.c)).*sqrt((2.*cons.Na.*cons.k.*T.*log(2))./(cons.M));
end


% ------------------------------
% Lorentzian pressure-broaddened
% at a certain p-T is 
%  ------------------------------

cont = 0;
for i=O2B_end:O2B_ini
    cont = cont+1;
    aux = (cons.T_ref./T).^(lines.temperatureDependence(i)); 
    gamma_p_T_O2B(cont,:) = aux.*((lines.airBroadenedWidth(i).*(p-p_self))+ (lines.selfBroadenedWidth(i).*p_self));
end

cont = 0;
for i=O2A_end:O2A_ini
    cont = cont+1;
    aux = (cons.T_ref./T).^(lines.temperatureDependence(i)); 
    gamma_p_T_O2A(cont,:) = aux.*((lines.airBroadenedWidth(i).*(p-p_self))+ (lines.selfBroadenedWidth(i).*p_self));
end




% ------------------------------
% Pressure shift-correction  
% of line position
%  ------------------------------
cont = 0;
for i=O2B_end:O2B_ini
    cont = cont+1;
    v_ij_O2B(cont,:) = lines.transitionWavenumber(i) + lines.pressureShift(i).*p;
end
cont = 0;
for i=O2A_end:O2A_ini
    cont = cont+1;
    v_ij_O2A(cont,:) = lines.transitionWavenumber(i) + lines.pressureShift(i).*p;
end
% ------------------------------
% Lower atmosphere: pressure broadening 
% of spectral lines dominates (assuming a Lorentz profile)
%  ------------------------------
cont = 0;
for i=O2B_end:O2B_ini
    cont = cont+1;
    f_L_O2B(cont,:,:) =  gamma_p_T_O2B(cont,:)./(pi.*(gamma_p_T_O2B(cont,:).^2 + (v_hr_O2B(cont,:)'-v_ij_O2B(cont,:)).^2)); % DIM: [line_i,HR,T]
end
cont = 0;
for i=O2A_end:O2A_ini
    cont = cont+1;
    f_L_O2A(cont,:,:) =  gamma_p_T_O2A(cont,:)./(pi.*(gamma_p_T_O2A(cont,:).^2 + (v_hr_O2A(cont,:)'-v_ij_O2A(cont,:)).^2)); % DIM: [line_i,HR,T]
end



% ------------------------------
% Upper atmosphere: Doppler boradening dominates
% of spectral lines dominates (assuming a Gaussian profile)
%  ------------------------------
cont = 0;
for i=O2B_end:O2B_ini
    cont = cont+1;
    f_G_O2B(cont,:,:) =  sqrt(log(2)./(pi*alpha_D_O2B(cont,:).^2)).*exp(-((v_hr_O2B(cont,:)'-v_ij_O2B(cont,:)).^2.*log(2))./(alpha_D_O2B(cont,:).^2)); % DIM: [line_i,HR,T]
end
cont = 0;
for i=O2A_end:O2A_ini
    cont = cont+1;
    f_G_O2A(cont,:,:) =  sqrt(log(2)./(pi*alpha_D_O2A(cont,:).^2)).*exp(-((v_hr_O2A(cont,:)'-v_ij_O2A(cont,:)).^2.*log(2))./(alpha_D_O2A(cont,:).^2)); % DIM: [line_i,HR,T]
end


% ------------------------------
% TEMPERATURE dependency on the 
% line intensity
%  ------------------------------
cont = 0;
for i=O2B_end:O2B_ini
    cont = cont+1;
    aux_1 = lines.lineIntensity(i);
    aux_2 = cons.Q./Q_T;
    aux_3 = exp(-cons.c2.*lines.lowerStateEnergy(i)./T)./exp(-cons.c2.*lines.lowerStateEnergy(i)./cons.T_ref);
    aux_4 = (double(1)-exp(-cons.c2.*v_ij_O2B(cont,:)'./T))./(double(1)-exp(-cons.c2.*v_ij_O2B(cont,:)'./cons.T_ref));
    S_ij_T_O2B(cont,:) = aux_1.*aux_2.*aux_3.*aux_4;
end

cont = 0;
for i=O2A_end:O2A_ini
    cont = cont+1;
    aux_1 = lines.lineIntensity(i);
    aux_2 = cons.Q./Q_T;
    aux_3 = exp(-cons.c2.*lines.lowerStateEnergy(i)./T)./exp(-cons.c2.*lines.lowerStateEnergy(i)./cons.T_ref);
    aux_4 = (double(1)-exp(-cons.c2.*v_ij_O2A(cont,:)'./T))./(double(1)-exp(-cons.c2.*v_ij_O2A(cont,:)'./cons.T_ref));
    S_ij_T_O2A(cont,:) = aux_1.*aux_2.*aux_3.*aux_4;
end


% ------------------------------
% TEMPERATURE dependency on the 
% line intensity PLOTTING 
%  ------------------------------
if flag_plot
    figure;
    cont = 0;
    cmap = redblue(length(T));
    colormap(cmap)
    
    for i=O2B_end:O2B_ini
        plot([nm(i) nm(i)],[lines.lineIntensity(i),0],'k') % Area
        hold on
    end
    
    for i=O2B_end:O2B_ini
        cont = cont+1;
        for tt= 1:length(T)
            plot(nm(i),S_ij_T_O2B(cont,tt),'o','Markersize',20,'color',cmap(tt,:)) % Maximun for a Lorentzian function (1/pi*gamma_HI)
            hold on
        end
    end
    %colormap(cmap(1:length(T),:))
   delta_T_med = (max(T)-min(T))./2;
   colorbar('Ticks',[0,0.5,1],...
        'TickLabels',{[num2str(min(T)-273),'^oC'],[num2str(min(T)-273+delta_T_med),'^oC'],[num2str(max(T)-273),'^oC']})

    set(gca,'Fontsize',18)
    title({'\Delta S_{ij} [cm^{-1}/molecule·cm^{-2}]'})
    xlabel('Wavelength [nm]')
    ylabel('S_{ij}=f(T)')
    set(gca,'TickLength', [0.03, 0.03])
    set(gca,'XMinorTick','on','YMinorTick','on')
    set(gcf,'PaperPositionMode','auto')
    set(gcf, 'Position', [300 500 800 400])
    
    if flag_save
        saveas(gcf, [figdir,filesep,'Delta_Sij_B_O16O',num2str(flags.iso),'.fig'],'fig');
        saveas(gcf, [figdir,filesep,'Delta_Sij_B_O16O',num2str(flags.iso),'.eps'],'epsc');
        saveas(gcf, [figdir,filesep,'Delta_Sij_B_O16O',num2str(flags.iso),'.png'],'png');
    end

       figure;
    cont = 0;
   
    
    for i=O2A_end:O2A_ini
        plot([nm(i) nm(i)],[lines.lineIntensity(i),0],'k') % Area
        hold on
    end
    
    for i=O2A_end:O2A_ini
        cont = cont+1;
        for tt= 1:length(T)
            plot(nm(i),S_ij_T_O2A(cont,tt),'o','Markersize',20,'color',cmap(tt,:)) % Maximun for a Lorentzian function (1/pi*gamma_HI)
            hold on
        end
    end
    colormap(cmap(1:length(T),:))
    colorbar('Ticks',[0,0.5,1],...
        'TickLabels',{[num2str(min(T)-273),'^oC'],[num2str(min(T)-273+delta_T_med),'^oC'],[num2str(max(T)-273),'^oC']})
    set(gca,'Fontsize',18)
    title({'\Delta S_{ij} [cm^{-1}/molecule·cm^{-2}]'})
    xlabel('Wavelength [nm]')
    ylabel('S_{ij}=f(T)')
    set(gca,'TickLength', [0.03, 0.03])
    set(gca,'XMinorTick','on','YMinorTick','on')
    set(gcf,'PaperPositionMode','auto')
    set(gcf, 'Position', [300 500 800 400])
    
    if flag_save
        saveas(gcf, [figdir,filesep,'Delta_Sij_A_O16O',num2str(flags.iso),'.fig'],'fig');
        saveas(gcf, [figdir,filesep,'Delta_Sij_A_O16O',num2str(flags.iso),'.eps'],'epsc');
        saveas(gcf, [figdir,filesep,'Delta_Sij_A_O16O',num2str(flags.iso),'.png'],'png');
    end
   
end


% ------------------------------
% COMPUTE the TRANSMITTANCE from the 
% Sij and fL (I can use the Lorentzian)
% bcs I am interested on the surface level
%  ------------------------------
% k_ij Monocromathic absorption coefficient f(T,P,\lambda)
% K_ij_O2A = S_ij_T_O2A.*f_L_O2A;
cont = 0;
for i=O2A_end:O2A_ini
    cont = cont+1;
    for tt= 1:length(T)
        K_ij_O2A(cont,:,tt) = S_ij_T_O2A(cont,tt).*f_L_O2A(cont,:,tt);
    end
end

cont = 0;
for i=O2B_end:O2B_ini
    cont = cont+1;
    for tt= 1:length(T)
        K_ij_O2B(cont,:,tt) = S_ij_T_O2B(cont,tt).*f_L_O2B(cont,:,tt);
    end
end


% Transmittance for each step of the INPUT file
% The volume density is [X]=[molecules/cm^3]  and u = [molecules/cm^2]
% If I assume that [X] is constant along the step then 
% u = [X]*L. L is the step in cm
% same value is constant in a step (e.g. 1 meter); I have to compute the 
% molecules in that step

if length(z_mes)==1 % If there is only one measure available (not a profile)
    % If there is only one measurement it is assumed to be the length of
    % the optical path on a nadir looking configuration
    z_step = z_mes.*100; % originally in [m] and converted to [cm]
    % O2rho is in [molecules/cm^3]
    % O2_den is in [molecules/cm^2]
    O2_den = O2rho.*z_step;
        
else % If there is a profile
    z_step = diff(z_mes).*100; % originally in [m] and converted to [cm]
    % O2rho is in [molecules/cm^3]
    % O2_den is in [molecules/cm^2]
    O2_den = (O2rho(1:end-1)+(diff(O2rho)./2)).*z_step;
end


% ------------------------------
% Compute the output spectral resolution
% ------------------------------

% MODTRAN resolution = 0.01cm-1
if strcmp(flags.resolution,'cm')
    wv_number_ini    = 1e7./nm(O2B_ini);
    wv_number_end    = 1e7./nm(O2B_end);
    wv_number_vector = wv_number_end:flags.res_num:wv_number_ini;
    nm_hr_O2B_vec    = 1e7./wv_number_vector;
    wl_O2B           = nm_hr_O2B_vec;
    
% resolution 0.01nm
elseif strcmp(flags.resolution,'nm')
    nm_hr_O2B_vec = sort([nm(O2B_ini):flags.res_num:nm(O2B_end)],'descend');
    wl_O2B        = nm_hr_O2B_vec(1:end-1)+(diff(nm_hr_O2B_vec)./2);
else 
    disp('No spectral resolution indicated')
end

% MODTRAN resolution = 0.01cm-1
if strcmp(flags.resolution,'cm')
    wv_number_ini    = 1e7./nm(O2A_ini);
    wv_number_end    = 1e7./nm(O2A_end);
    wv_number_vector = wv_number_end:flags.res_num:wv_number_ini;
    nm_hr_O2A_vec    = 1e7./wv_number_vector;
    wl_O2A           = nm_hr_O2A_vec;
    % resolution 0.01nm
elseif strcmp(flags.resolution,'nm')
    nm_hr_O2A_vec = sort([nm(O2A_ini):flags.res_num:nm(O2A_end)],'descend');
    wl_O2A        = nm_hr_O2A_vec(1:end-1)+(diff(nm_hr_O2A_vec)./2);
else
    disp('No spectral resolution indicated')
end


% ------------------------------
% O2 RAYLEIGH OPTICAL depth
% For the moment it has been computed the 
% O2 optical depth associated with the O2 absorption TAU_GAS(O_2^16,O16O17,O16O18) 
% However, the optical depth associated with the Rayleigh scattering of the
% O2 should be also computed. In general, TAU = TAU_GAS + TAU_RAY + TAU_AER
% ------------------------------


% This expression uses the wavelength in micrometers wl_O2B*1e-3 and wl_O2A*1e-3
wl_O2B = wl_O2B./1e3; wl_O2A = wl_O2A./1e3;

scatte_cross_sec_2B= 1e-28.*(1.0455996-341.29061.*wl_O2B.^(-2)-0.90230850.*wl_O2B.^2)./(1+0.0027059889.*wl_O2B.^(-2)-85.968563.*wl_O2B.^(2)); %[cm2] % Eq. 29 from
scatte_cross_sec_2A= 1e-28.*(1.0455996-341.29061.*wl_O2A.^(-2)-0.90230850.*wl_O2A.^2)./(1+0.0027059889.*wl_O2A.^(-2)-85.968563.*wl_O2A.^(2)); %[cm2] % Eq. 29 from
% approximation used in: On Rayleigh Optical Depth Calculations
% BARRY A. BODHAINE
wl_O2B = wl_O2B.*1e3; wl_O2A = wl_O2A.*1e3;


tau_ij_ray_O2B        = scatte_cross_sec_2B.*p.*cons.Na./(cons.ma.*cons.g); % Rayleigh contribution
tau_ij_ray_O2A        = scatte_cross_sec_2A.*p.*cons.Na./(cons.ma.*cons.g); % Rayleigh contribution

% Units: p is in [atm], scatte_cross_sec_ [cm2], cons.Na [mol-1], cons.g[m·s-2], cons.ma [g·mol-1]
% So, to make the tau_ij_ray_O2A adimentional we need to use the factor to
% convert 1 atm to 101325 Pa, and 10 from the m2->cm2 and g->kg
tau_ij_ray_O2B        = tau_ij_ray_O2B.*101325./10;
tau_ij_ray_O2A        = tau_ij_ray_O2A.*101325./10;


tran_ray_O2B          = exp(-tau_ij_ray_O2B./cosd(zenith));
tran_ray_O2A          = exp(-tau_ij_ray_O2A./cosd(zenith));

cont = 0;
for i=O2B_end:O2B_ini
    cont = cont+1;
    for tt= 1:length(T)
        tau_ij_O2B            = O2_den(tt).*K_ij_O2B(:,:,tt); % Absorption contribution   
        tran_O2B(cont,:,tt)   = exp(-tau_ij_O2B(cont,:)./cosd(zenith));
       
    end
end


cont = 0;
for i=O2A_end:O2A_ini
    cont = cont+1;
    for tt= 1:length(T)
        tau_ij_O2A            = O2_den(tt).*K_ij_O2A(:,:,tt); % Absorption contribution
        tran_O2A(cont,:,tt)   = exp(-tau_ij_O2A(cont,:)./cosd(zenith));
      
    end
end

Total_tran_A = prod(tran_O2A,3);
Total_tran_B = prod(tran_O2B,3);


Total_tran_A_ray = prod(tran_ray_O2A,1);
Total_tran_B_ray = prod(tran_ray_O2B,1);



% Grouping and convolution 
FWHM       = abs(1e7/((1e7/690)+0.01)-690);
step_scale = 50;
for i=1:length(nm_hr_O2B_vec)
    
    % select the wavelength range    
    if i==length(nm_hr_O2B_vec) % last wv
        delta = abs(nm_hr_O2B_vec(i)-nm_hr_O2B_vec(i-1))./2;
    else
        delta = abs(nm_hr_O2B_vec(i)-nm_hr_O2B_vec(i+1))./2;
    end
    
    for j=1:size(nm_hr_O2B,1)
        [row,col]      = find(nm_hr_O2B(j,:)<=(nm_hr_O2B_vec(i)+delta) & nm_hr_O2B(j,:)>=(nm_hr_O2B_vec(i)-delta));   
        % Obtain the information of each wavelength interval and for each transition line
        if ~isempty(col) && numel(col)~=1
            Total_tran_B_inter(j,:) = interp1(nm_hr_O2B(j,col),Total_tran_B(j,col),[nm_hr_O2B_vec(i)-step_scale.*delta:1e-4:nm_hr_O2B_vec(i)+step_scale.*delta],'linear','extrap');   
        elseif ~isempty(col) && numel(col)==1
            Total_tran_B_inter(j,:) = repmat(Total_tran_B(j,col),1,length([nm_hr_O2B_vec(i)-step_scale.*delta:1e-4:nm_hr_O2B_vec(i)+step_scale.*delta]));
        else 
            Total_tran_B_inter(j,:) = ones(1,length([nm_hr_O2B_vec(i)-step_scale.*delta:1e-4:nm_hr_O2B_vec(i)+step_scale.*delta]));
        end 
    end
    
    % All the transition lines contribution 
    TRANS_row = prod(Total_tran_B_inter,1);

    % Define the ISRF
    sig = FWHM/2.3548;
    C   = nm_hr_O2B_vec(i);
    % Calculate ISRF:
    wvl = [nm_hr_O2B_vec(i)-step_scale.*delta:1e-4:nm_hr_O2B_vec(i)+step_scale.*delta];
    isrf = exp(-(wvl - C).^2/(2*sig^2));    
   
    % sig: std deviation
    x   = wvl;
    mu  = C; % mu : mean
    amp = 1; % amp: negative or postive
    vo  = 0; % vo: vertival offset from the baseline
    
    gaus = @(x,mu,sig,amp,vo)amp*exp(-(((x-mu).^2)/(2*sig.^2)))+vo;
    isrf = gaus(x,mu,sig,amp,vo);
    
    if 1==0
        figure
        plot(wvl,isrf)
    end
    
    % convolve the signal
    isrf_norm = isrf./trapz(wvl,isrf);
    TRAN_O2B(i) = trapz(wvl,(TRANS_row.*isrf_norm)',1);
   
    clear Total_tran_B_inter TRANS_row     
end


% -----
% O2-A 
% -----

% Grouping and convolution 
FWHM = abs(1e7/((1e7/765)+0.01)-765);
for i=1:length(nm_hr_O2A_vec)
    
    % select the wavelength range
    if i==length(nm_hr_O2A_vec) % last wv
        delta = abs(nm_hr_O2A_vec(i)-nm_hr_O2A_vec(i-1))./2;
    else
        delta = abs(nm_hr_O2A_vec(i)-nm_hr_O2A_vec(i+1))./2;
    end
    
    for j=1:size(nm_hr_O2A,1)
        [row,col]      = find(nm_hr_O2A(j,:)<=(nm_hr_O2A_vec(i)+delta) & nm_hr_O2A(j,:)>=(nm_hr_O2A_vec(i)-delta));   
        % Obtain the information of each wavelength interval and for each transition line
        if ~isempty(col) && numel(col)~=1
            Total_tran_A_inter(j,:) = interp1(nm_hr_O2A(j,col),Total_tran_A(j,col),[nm_hr_O2A_vec(i)-step_scale.*delta:1e-4:nm_hr_O2A_vec(i)+step_scale.*delta],'linear','extrap');   
        elseif ~isempty(col) && numel(col)==1
            Total_tran_A_inter(j,:) = repmat(Total_tran_A(j,col),1,length([nm_hr_O2A_vec(i)-step_scale.*delta:1e-4:nm_hr_O2A_vec(i)+step_scale.*delta]));
        else 
            Total_tran_A_inter(j,:) = ones(1,length([nm_hr_O2A_vec(i)-step_scale.*delta:1e-4:nm_hr_O2A_vec(i)+step_scale.*delta]));
        end 
    end
    
    % All the transition lines contribution 
    TRANS_row = prod(Total_tran_A_inter,1);

    % Define the ISRF
    sig = FWHM/2.3548;
    C   = nm_hr_O2A_vec(i);
    % Calculate ISRF:
    wvl = [nm_hr_O2A_vec(i)-step_scale.*delta:1e-4:nm_hr_O2A_vec(i)+step_scale.*delta];
    isrf = exp(-(wvl - C).^2/(2*sig^2));    
   
    % sig: std deviation
    x   = wvl;
    mu  = C; % mu : mean
    amp = 1; % amp: negative or postive
    vo  = 0; % vo: vertival offset from the baseline
    
    gaus = @(x,mu,sig,amp,vo)amp*exp(-(((x-mu).^2)/(2*sig.^2)))+vo;
    isrf = gaus(x,mu,sig,amp,vo);
    
    if 1==0
        figure
        plot(wvl,isrf)
    end
    
    % convolve the signal
    isrf_norm = isrf./trapz(wvl,isrf);
    TRAN_O2A(i) = trapz(wvl,(TRANS_row.*isrf_norm)',1);
   
    clear Total_tran_A_inter TRANS_row     
end






if flag_plot
    figure; 
    plot(wl_O2B,TRAN_O2B,'r-')
    hold on
    plot(nm_hr_O2B,Total_tran_B,'k.')
    plot(wl_O2B,TRAN_O2B,'r-','Linewidth',2)
    legend('SR=0.01mn & SSI = 0.01 nm','Original HITRAN resolution')
    xlabel('Wavelength [nm]')
    set(gca,'Fontsize',18)
    ylabel('Transmittance [-]')
    set(gca,'TickLength', [0.03, 0.03])
    set(gca,'XMinorTick','on','YMinorTick','on')
    
    title('O_2-B region')
    if flag_save
        saveas(gcf, [figdir,filesep,'Total_tran_B_O16O',num2str(flags.iso),'.fig'],'fig');
        saveas(gcf, [figdir,filesep,'Total_tran_B_O16O',num2str(flags.iso),'.eps'],'epsc');
        saveas(gcf, [figdir,filesep,'Total_tran_B_O16O',num2str(flags.iso),'.png'],'png');
    end
        
    figure
    plot(wl_O2A,TRAN_O2A,'r-')
    hold on
    plot(nm_hr_O2A,Total_tran_A,'k.')
    plot(wl_O2A,TRAN_O2A,'r-','Linewidth',2)
    legend('SR=0.01mn & SSI = 0.01 nm','Original HITRAN resolution')
    set(gca,'Fontsize',18)
    xlabel('Wavelength [nm]')
    ylabel('Transmittance [-]')
    set(gca,'TickLength', [0.03, 0.03])
    set(gca,'XMinorTick','on','YMinorTick','on')

    title('O_2-A region')
    
    if flag_save
        saveas(gcf, [figdir,filesep,'Total_tran_A_O16O',num2str(flags.iso),'.fig'],'fig');
        saveas(gcf, [figdir,filesep,'Total_tran_A_O16O',num2str(flags.iso),'.eps'],'epsc');
        saveas(gcf, [figdir,filesep,'Total_tran_A_O16O',num2str(flags.iso),'.png'],'png');
    end
end



 
end 
