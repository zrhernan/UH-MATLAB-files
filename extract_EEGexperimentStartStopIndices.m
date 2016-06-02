clc, clear all, close all
tic  % start computation time

% for calling helper functions
addpath(genpath('\\bmi-nas-01\Contreras-UH\Infantdata\Infantdata\code\Zachs_Infant_decoding_files'))

%% Annotation for Infant Data
files = dir('\\bmi-nas-01\Contreras-UH\Infantdata\Infantdata\Data\');   % assume starting from current directory
filenames = {files.name};
subdirs = filenames([files.isdir]);
expr = '-\d{4}$'; % contains the year at the end of the folder name
cnt = 1;
for s = 1:length(subdirs)
    if ~isempty(regexp(subdirs{s},expr, 'once'))
        InfantDataAnnotList{cnt} = subdirs{s};
        hypenIDX = regexp(subdirs{s},'-');
        InfantID{cnt} = subdirs{s}(1:hypenIDX(1)-1);
%         disp(['Infant folder recognized: ',InfantDataAnnotList{cnt}])
        cnt = cnt + 1;
    end
end


% infant = menu(['Which infant would you',10,' like to analyze?'],InfantID{:});
% close
%%
for infant = 1:length(InfantID);
%% Initializing and Assigning directory paths....
serverPath1 = ['\\bmi-nas-01\Contreras-UH\Infantdata\Infantdata\Data\',...
    InfantDataAnnotList{infant},'\EEG'];
% serverPath2 = ['\\bmi-nas-01\Contreras-UH\Infantdata\Infantdata\code\',...
%     'Zachs_Infant_decoding_files\Data\',InfantDataAnnotList{infant}];
% desktopPath = ['C:\Users\zrhernan\Infant_decoding_files\Data\',...
%     InfantDataAnnotList{infant}];
%
%% Check if Experiment Start-Stop Indices have been saved already
fullFileName1 = [serverPath1,'\ExperStartStopIndices.mat'];
fullFileName2 = [serverPath1,'\EEGTriggerTimes.mat'];
if exist(fullFileName1, 'file') || exist(fullFileName2, 'file')
  % File already exists.
  if exist(fullFileName1, 'file') 
    warningMessage = sprintf('Warning: ''ExperStartStopIndices'' file already exists:\n%s', fullFileName1);
    disp(warningMessage)
  elseif exist(fullFileName2, 'file')
    warningMessage = sprintf('Warning: ''EEGTriggerTimes'' file already exists:\n%s', fullFileName1);
    disp(warningMessage)
  end
  disp('Skipping to next infant data set')
  continue
end

%% Check if EEG '.vhdr' file exists as well
%
fullFileName = [serverPath1,'\',InfantDataAnnotList{infant},'.vhdr'];
if ~exist(fullFileName,'file')
    % File does not exist.
    warningMessage = sprintf('Warning: EEG VHDR file does not exist:\n%s', fullFileName);
    disp(warningMessage)
    disp('Skipping to next infant data set')
    continue
end

%% Opening the VHDR file
EEG = pop_loadbv(serverPath1,[InfantDataAnnotList{infant},'.vhdr']);
disp('____Brain Vision Recorder EEG file imported')
disp('--------------------------------------------------')
timecompute(toc)
disp('--------------------------------------------------')

%% Opening Event Markers (triggers) and plotting them as a bar graph
eventTriggers = zeros(1,length(EEG.event));
for trig = 1:length(EEG.event)
    eventTriggers(trig) = EEG.event(trig).latency;
end
bar(eventTriggers)

%% Input the Experiment Start and Stop Times
ExperStartInd = input(['Based on the bar plot, ', 10,...
    'please input the index corresponding to the start of the experiment: ']);
ExperEndInd = input(['Based on the bar plot, ', 10,...
    'please input the index corresponding to the end of the experiment: ']);
%}

%% Saving the resulting variables as a MAT file
    cd(serverPath1) % save to Infantdata path
    
    filename2save1 = 'ExperStartStopIndices';
    save([filename2save1,'.mat'],'ExperStartInd','ExperEndInd');
%     load([serverPath1,'\',filename2save1,'.mat'])  % Loading Experiment Start-Stop Times (did this cuz I forgot to save in a text file)
    fileID = fopen([filename2save1,'.txt'],'w');
    fprintf(fileID,'%d\n',ExperStartInd);
    fprintf(fileID,'%d\n',ExperEndInd);
    fclose(fileID);
    
    filename2save2 = 'EEGTriggerTimes';
    save([filename2save2,'.mat'],'eventTriggers');
    fileID = fopen([filename2save2,'.txt'],'w');
    fprintf(fileID,'%d\n',eventTriggers);
    fclose(fileID);
    
    disp('____experiment start-stop indices for EEG data saved to lab server')
%     cd(desktopPath) % save to path on Zach's desktop
%     save(filename2save,'ExperStartInd','ExperEndInd');
%     disp('____EEG data structure saved to local Desktop')

disp('--------------------------------------------------')
timecompute(toc)
disp('--------------------------------------------------')
close all
end  % repeat for each infantID