%% For generating a class label vector from start and stop time segments of each class
clc, clear all, close all
tic;
% for calling helper functions
addpath('\\172.27.216.40\Contreras-UH\Infantdata\Infantdata\code\Zachs_Infant_decoding_files') 
%% Annotation for Infant Data
InfantDataAnnotList = {'N09-12-04-2013','J20-08-19-2013','B06-10-30-2013',...
    'GR09-07-12-2014','A06-09-28-2013','LW10-06-19-2014','AR16-07-14-2014',...
    'RB23-12-04-2014','A18-08-15-2013'};
% LastColumnIndex_7C_J20 = 162;  % for the seven-class dataset of subject J20
LastColumnIndex = [104,144,68,114,49,50,153,238,83]; %for multiple class dataset
% LastColumnIndex = [39,46,19]; %for binary class dataset
infant = menu(['Which infant would you',10,' like to analyze?'],'N09',...
    'J20','B06','GR09','A06','LW10','AR16','RB23','A18');
close

%% Path directory initializations
serverPath1 = ['\\172.27.216.40\Contreras-UH\Infantdata\Infantdata\Data\',...
    InfantDataAnnotList{infant}];
serverPath2 = ['\\172.27.216.40\Contreras-UH\Infantdata\Infantdata\code\Zachs_Infant_decoding_files\Data\',...
    InfantDataAnnotList{infant}];
desktopPath = ['C:\Users\zrhernan\Infant_decoding_files\Data\',...
    InfantDataAnnotList{infant}];

%% Importation of Class Start and Stop Times
% filename = [desktopPath,'\Observe&Imitate_start-stop times.txt'];
% filename = [desktopPath,'\7-class start-stop times.txt'];
filename = [desktopPath,'\class start-stop times.txt'];
TASKSEGMENTS = importTaskTrialInfo(filename, 2, LastColumnIndex(infant));   %_7C_J20);
tasklabels = capitalize(unique(TASKSEGMENTS.Task));     % create a list of all task labels

%% Importation of EEG Data
load([desktopPath,'\EEGfiles.mat'],'timesampleStart','timesampleEnd');  % load structure of EEG attributes

%% Shifting Data to match timing of video recording
% length of video synced recording in EEG {under Fs = 100 Hz}
EEGtimelength = round((timesampleEnd - timesampleStart)/10); 

%% Converting to 100 Hz sampling frequency
TASKSEGMENTS.StartTime = round((TASKSEGMENTS.StartTime)/10);
TASKSEGMENTS.StopTime = round((TASKSEGMENTS.StopTime)/10);

%% Generate the class label vector for further classification
classlabel = zeros(EEGtimelength,1);
for triali = 1:size(TASKSEGMENTS.TaskLabel,1)
    classlabel(TASKSEGMENTS.StartTime(triali):TASKSEGMENTS.StopTime(triali)) = TASKSEGMENTS.TaskLabel(task);
end

%% Plot the labeling vector vs. time samples
classtargetvector_plot(classlabel, tasklabels);

%% Save the classlabel file
cd(serverPath1) % save to Infantdata path where 'Data' is located
save('classlabel.mat', 'classlabel');
disp('____classlabel vector saved to folder "Data" lab server')
disp('--------------------------------------------------')
timecompute(toc)
disp('--------------------------------------------------')
cd(serverPath2) % save to Infantdata path where 'Zachs_Infant_decoding_files' are located
save('classlabel.mat', 'classlabel');
disp('____classlabel vector saved to lab server folder "Zachs_Infant_decoding_files"')
disp('--------------------------------------------------')
timecompute(toc)
disp('--------------------------------------------------')
cd(desktopPath) % save to path on Zach's desktop
save('classlabel.mat', 'classlabel');
disp('____classlabel vector saved to local Desktop')
disp('--------------------------------------------------')
timecompute(toc)
disp('--------------------------------------------------')