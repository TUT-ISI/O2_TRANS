function [num_samples,AUX,AUX_GEO,path] = read_inputs()
% READ_INPUTS This function ask for the input (*.txt) file where data of measured
% pressure, Temperature and surface elevation is contained.

% INPUT: 
%     INPUTfile.txt must contain at least three columns being 
%       1st column Z [m]
%       2nd column p [mbar]
%       3rd column T [K]
%     The each row corresponds to a different measurement taken at a different elevation
%     If more than one sample wants to be evaluated (e.g. measurements
%     different days) attach them as extra columns always (Z,p,T) comma
%     separated.

% OUTPUT: 
%     num_samples: This files reads the number of triplet columns (Z.p,T)
%     available in the INPUT file. According to the num_samples we generate
%     a loop
% 
% Author: Neus Sabater
% Version v.0
% Data: April/2020
% e-mail: neus.sabater@fmi.fi
% ----------------------------------------------------------------------------------

if ~ispc %To force showing a dialogue in MAC and UNIX os
     questdlg('Please, pick the input file containing (Z[m], p[mbar], T[K])', ...
        'INPUT file',...
        'OK','OK'); 
end



[file, path, ~] = uigetfile('*.txt', 'Please, pick the input file containing (Z[m],p[mbar],T[K])');


AUX = importdata([path,file]);


num_samples = size(AUX,2)/3;

if mod(abs(num_samples),1)~=0
    disp(['Wrong number of columns in the INPUT file: ',file])
    return 
end 



if ~ispc %To force showing a dialogue in MAC and UNIX os
     questdlg('Please, pick the geometry file containing the solar or viewing zenith angle in degrees for each profile', ...
        'INPUT file',...
        'OK','OK'); 
end



[file, path, ~] = uigetfile('*.txt', 'pick the geometry file that contain the solar or viewing zenith angle in degrees for each profile');


AUX_GEO = importdata([path,file]);


num_samples_GEO = size(AUX_GEO,2);

if num_samples_GEO~=num_samples
    disp('Different number of columns in the INPUT and in the GEOMETRY files')
    return 
end 






end
