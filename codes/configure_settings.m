function flags= configure_settings()
% CONFIGURE_SETTINGS This function ask for settings configuration.
% Particularly, it ask the user if intermedaited results needs to be
% plotted; figures needs to be saved and in that case the path where to
% store the figures

% SETTINGS: 

%     flags.plotting = 0; [0=no plotting] [1=plotting]
%     flags.saving = 0;   [0=no saving] [1=saving]
%     flags.dir = directory; Only required if flags.saving =1

% Author: Neus Sabater
% Version v.0
% Data: April/2020
% e-mail: neus.sabater@fmi.fi
% ----------------------------------------------------------------------------------


if ~ispc
     questdlg('Pick the directory of the AUX folder', ...
        'AUX folder', ...
        'OK','OK'); 
end

flags.AUX      = uigetdir('Pick the directory of the AUX folder');


answer = questdlg('Do you want to plot intermediate results?', ...
    'Plotting options', ...
    'YES','NO','NO');
% Handle response
switch answer
    case 'YES'
        flags.plotting =1;
        
    case 'NO'
        flags.plotting =0;
end


% flags.plotting = input('Do you want to plot intermediate results? YES = [1]; NO=[0]');


if  flags.plotting ==1
    answer = questdlg('Do you want to save plots (*.png,*.eps,*.fig)?', ...
        'Plotting options', ...
        'YES','NO','NO');
    % Handle response
    switch answer
        case 'YES'
            flags.saving = 1;
            
        case 'NO'
            flags.saving = 0;
    end
    if flags.saving==1
        if ~ispc %To force showing a dialogue in MAC and UNIX os
            questdlg('Please, Pick the directory to save the figures', ...
                'Saving directory', ...
                'OK','OK');
        end
        flags.dir = uigetdir;
    else
        flags.dir =[];
    end
else
    flags.saving =0;
    flags.dir =[];
end

% Ask about the resolution for the output numbers
warning off
answer = questdlg('Would you like to define the spectral resolution and interval in...?', 'Spectral resolution',...
	'nm','cm','default');

prompt = {'Enter spectral resolution '};
definput = {'0.1'};
dlgtitle = 'Input';
dims     = [1 35];
% Handle response
switch answer
    case 'nm'
        flags.resolution = 'nm';
        aux = inputdlg(prompt,dlgtitle,dims,definput);
        flags.res_num = str2num(aux{1});
    case 'cm'
        flags.resolution = 'cm';
        aux = inputdlg(prompt,dlgtitle,dims,definput);
        flags.res_num = str2num(aux{1});
end


end 
