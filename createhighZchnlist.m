function [ highZchnlist, Zchnlist ] = createhighZchnlist( filepaths )
%CREATEHIGHZCHNLIST imports impedance values using the 
%   IMPORTEEGIMPEDANCES.M file and consolidates a list of channels with
%   high impedances before and after the experiment
%   INPUT: filepaths        = file paths to the 'Brain Products actiCAP Impedance 
%                             File Version 1.0' text file to import (at
%                             start and end of experiment)
%   OUTPUT: highZchnlist    = list of channels with high impedances (Z)
%                             recorded before and after the experiment
%   created by: Zach Hernandez, University of Houston, 2015

%% Check if Experiment Start-End Impedance Files exist
if ~exist(filepaths{1}, 'file') && ~exist(filepaths{2}, 'file')
    % File does not exist.
    warningMessage = sprintf('Warning: Both impedance files do not exist');
    disp(warningMessage)
    highZchnlist = [];
    Zchnlist = {nan(64,1),nan(64,1)};
    return  % break out of function
elseif ~exist(filepaths{1}, 'file') || ~exist(filepaths{2}, 'file')
    if ~exist(filepaths{1}, 'file')
        % File does not exist.
        warningMessage = sprintf('Warning: file does not exist:\n%s', filepaths{1});
        disp(warningMessage)
        impedancefilePaths = {'';filepaths{2}};
    elseif ~exist(filepaths{2}, 'file')
        % File does not exist.
        warningMessage = sprintf('Warning: file does not exist:\n%s', filepaths{2});
        disp(warningMessage)
        impedancefilePaths = {filepaths{1};''};
    end
else
    impedancefilePaths = filepaths;
end


%% IMPEDANCE FILE FROM START OF EXPERIMENT
if ~strcmp(impedancefilePaths{1},'') % check if impedance filename exists
% 1) import Brain Products impedance files    
    [EL_start, IV_start] = importEEGImpedances(impedancefilePaths{1}); % recorded at beginning of experiment
% 2) extract channels with high impedances(displayed as 'NaN's in list of
% impedances
    highZlistSTART = EL_start(isnan(IV_start) | IV_start > 60);
    IV_start(isnan(IV_start)) = 100;
else
    highZlistSTART = {};
    IV_start = nan(64,1);
end

%% IMPEDANCE FILE FROM END OF EXPERIMENT
if ~strcmp(impedancefilePaths{2},'') % check if impedance filename exists
% 1) import Brain Products impedance files
    [EL_end, IV_end] = importEEGImpedances(impedancefilePaths{2});     % recorded at end of experiment
% 2) extract channels with high impedances(displayed as 'NaN's in list of
% impedances    
    highZlistEND = EL_end(isnan(IV_end) | IV_start > 60);
    IV_end(isnan(IV_end)) = 100;
else
    highZlistEND = {};
    IV_end = nan(64,1);
end


% 3) combined both list channels together without any overlap
highZchnlist = transpose(unique([highZlistSTART; highZlistEND]));
Zchnlist = {IV_start, IV_end};

end

