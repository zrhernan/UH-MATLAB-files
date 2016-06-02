function [ InfantDataAnnotList, InfantID ] = defineInfantFolders( InfantDir, disp_flag )
%DEFINEINFANTFOLDERS generates List of infant data folders and the infant
%ID name for labeling and indexing purposes.
% Input:    InfantDir -->   folder where all infant data is located 
%                           (includes video, segmentation, EEG, and 
%                           kinematics files for each infant (subject))
%           disp_flag -->   decides whether or not to display the infant 
%                           subject folder being recognized by the function
%
% Output:   InfantDataAnnotList     --> list of folders for each infant
%                                       (1-by-size of recognized subjects)
%           InfantID                --> list of indexing labels for each
%                                       infant (same size as
%                                       'InfantDataAnnotList')
%   Created by: Zach Hernandez, University of Houston, 2016
%% Generate List of Infant Data Folders
if nargin < 1
    InfantDir = '\\172.27.216.40\Contreras-UH\Infantdata\Infantdata\Data\';
elseif nargin < 2
    disp_flag = 1;
end

files = dir(InfantDir);   % assume starting from current directory
filenames = {files.name};
subdirs = filenames([files.isdir]);
expr = '-\d{4}$'; % contains the year at the end of the folder name
cnt = 1;
for s = 1:length(subdirs)
    if ~isempty(regexp(subdirs{s},expr, 'once'))
        InfantDataAnnotList{cnt} = subdirs{s};
        hypenIDX = regexp(subdirs{s},'-');
        InfantID{cnt} = subdirs{s}(1:hypenIDX(1)-1);
        if disp_flag, disp(['Infant folder recognized: ',InfantDataAnnotList{cnt}]), end
        cnt = cnt + 1;
    end
end

end   %EOF