clc, clear all, close all
tic  % start computation time

% for calling helper functions
addpath(genpath('\\bmi-nas-01\Contreras-UH\Infantdata\Infantdata\code\Zachs_Infant_decoding_files'))

%% Generate List of Infant Data Folders
InfantDir = '\\172.27.216.40\Contreras-UH\Infantdata\Infantdata\Data\';

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
        disp(['Infant folder recognized: ',InfantDataAnnotList{cnt}])
        cnt = cnt + 1;
    end
end

%%
for infant = 1:61;
disp(['Opening EEG recording of subject ',InfantID{infant}])
  
%% Initializing and Assigning directory paths....
    infantfolder = InfantDataAnnotList{infant};
    serverPath1 = [InfantDir,infantfolder,'\Kinematics\'];
    cd(serverPath1)

%% Check if Experiment Start-Stop Indices have been saved already
fullFileName = ['ExperStartStopIndices.mat'];
% if exist(fullFileName, 'file')
%   % File does not exist.
%   warningMessage = sprintf('Warning: file already exists:\n%s', fullFileName);
%   disp(warningMessage)
%   disp('Skipping to next infant data set')
%   continue
% end

%% Sensor ID to Mounted Body Part List
sensorIDList_post2013 = {'SI-000719','SI-000775','SI-000722','SI-000773',...
    'SI-000738','SI-000708'}; % OPAL sensor IDs for experiments after 2013
sensorIDList_2013 = {'SI-000708','SI-000719','SI-000722','SI-000738',...
    'SI-000742','SI-000773'}; % OPAL sensor IDs for 2013 experiments
sensorBPList = { 'ExperimentersLeftArm','ExperimentersRightArm',...
    'InfantsLeftArm','InfantsRightArm','InfantsTrunk','InfantsForehead'};
sensorID = 3; % only time from one sensor needed (since timing is the same for all sensors) 

% Specify which sensor ID set to use (depends on the year experiment was conducted)
 if ~isempty(strfind(InfantDataAnnotList{infant},'2014')) || ~isempty(strfind(InfantDataAnnotList{infant},'2015')) || ~isempty(strfind(InfantDataAnnotList{infant},'2016'));
     sensorIDList = sensorIDList_post2013;
 elseif ~isempty(strfind(InfantDataAnnotList{infant},'2013'));
     sensorIDList = sensorIDList_2013;
 end
 disp('Reading IMU data for....')
 disp(['........Subject ' InfantID{infant}])
 disp(['........' sensorBPList{sensorID} ' | Sensor ID ' sensorIDList{sensorID}])

%% Set directory paths for related scripts
    addpath(genpath('\\bmi-nas-01\Contreras-UH\Lab software and hardware\quaternions'))

%% Read OPAL sensor data directly from HDF5 (.h5) file format. Refer Pg. 44 of OPAL sensor's Manual

filename = [InfantDataAnnotList{infant},'.h5'];

try
    vers = h5readatt(filename,'/','FileFormatVersion');
catch err
    try
        vers = h5readatt(filename,'/','File_Format_Version');
    catch err
        error('Couldn''t determine file format');
    end
end
if vers < 2
    error('This example only works with version 2 or later of the data file')
end

caseIdList = hdf5read(filename,'/CaseIdList');      % Difficult to convert to h5read

for SI = 1:size(caseIdList,1);
    groupName = ['/' caseIdList(SI).data]; 
    if strcmp(groupName,['/',sensorIDList{sensorID}]);
        break      % extract head sensor from .h5 file automatically
    else
        clear groupName  
    end
end
if ~exist('groupName','var');
    disp(['ERROR:',10,'   OPAL sensor ',sensorIDList{sensorID},' NOT found. Re-check available',10,' sensors for Subject ' InfantID{infant}]);
    return
end

Time_Path  = [groupName '/Time'];

Annotations = h5read(filename,'/Annotations');
Time_stamps = h5read(filename,Time_Path);

fs = h5readatt(filename, groupName, 'SampleRate');
fs = double(fs);    % Automatically reads the sampling frequency
Ts = 1/fs;
t_plot = (1:size(Time_stamps,1))/fs;

%% Acquiring Event Markers (triggers)from Annotation Time Stamps and plotting them as a bar graph
% if strcmp(Annotations.Annotation(27,1),'0')
%     Annotations.Time = Annotations.Time(2:end,:);
% end
if strcmp(filename,'NG16-11-14-2015.h5')        %just solving a confilct with NG16.... the trigger did not record with Opals. It did with EEG.
    a = Annotations.Time(end);
    b(1) = a+uint64(10); Annotations.Time(14) = uint64(b(1));
    for i = 1:10; b(i+1) = b(1)+uint64(1*60*10^6+10*i); Annotations.Time(i+14) = uint64(b(i+1)); end
end
if strcmp(filename,'VC06-06-26-2015.h5')        %just solving a confilct with VC06.... A random high trigger happened in the middle.
    Annotations.Time(9)=[];
end
FallEdgeIDX = find(Annotations.Annotation(27,:)=='+');
TrigEdgeOnset = double(Annotations.Time(FallEdgeIDX)-min(Time_stamps));
TrigEdgeOnset_128Hz = double(round(TrigEdgeOnset.*(128/1e+6)));  %re-format so that onset values are in 128 Hz samples
TrigEdgeOnset_1000Hz = double(round(TrigEdgeOnset.*(1000/1e+6)));  %re-format so that onset values are in 1000 Hz samples
bar(1:length(TrigEdgeOnset_1000Hz),TrigEdgeOnset_1000Hz)
set(gca,'xtick',1:length(TrigEdgeOnset_1000Hz))

 
%% Input the Experiment Start and Stop Times
load('ExperStartStopIndices.mat')
if ExperStartInd==0
    disp(['Start == 0 for ',InfantID{infant}])
end
if ExperEndInd==0
    disp(['End == 0 for ',InfantID{infant}])
end
% ExperStartInd = input(['Based on the bar plot, ', 10,...
% 'please input the index corresponding to the start of the experiment: ']);
% ExperEndInd = input(['Based on the bar plot, ', 10,...
% 'please input the index corresponding to the end of the experiment: ']); 

%% Saving the resulting variables as a MAT file
    cd(serverPath1) % save to Infantdata path
    
    filename2save1 = 'ExperStartStopIndices';
    save([filename2save1,'.mat'],'ExperStartInd','ExperEndInd');
%     load([serverPath1,'\',filename2save1,'.mat'])  % Loading Experiment Start-Stop Times (did this cuz I forgot to save in a text file)
    fileID = fopen([filename2save1,'.txt'],'w');
    fprintf(fileID,'%d\n',ExperStartInd);
    fprintf(fileID,'%d\n',ExperEndInd);
    fclose(fileID);
    
    filename2save2 = 'KINETriggerTimes';
    save([filename2save2,'.mat'],'TrigEdgeOnset_1000Hz');
    fileID = fopen([filename2save2,'.txt'],'w');
    fprintf(fileID,'%d\n',TrigEdgeOnset_1000Hz);
    fclose(fileID);
    
    disp('____experiment start-stop indices for kinematics data saved to lab server')


    disp('--------------------------------------------------')
    timecompute(toc)
    disp('--------------------------------------------------')
    close all
end  % repeat for each infantID


