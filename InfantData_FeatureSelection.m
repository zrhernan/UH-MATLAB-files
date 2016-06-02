clc, clear all, format compact, %close all
tic  % start computation time

%% for calling helper functions
addpath(genpath('\\172.27.216.40\Contreras-UH\Infantdata\Infantdata\code\Zachs_Infant_decoding_files'))

%% Annotation for Infant Data
InfantDataAnnotList = {
'A06-09-28-2013',...
'A07-08-12-2013',...
'A18-08-15-2013',...
'AD10-03-24-2015',...
'AG06-03-02-2015',...
'AH10-04-22-2015',...
'AR16-07-14-2014',...
'AS13-01-31-2015',...
'B06-10-30-2013',...
'BG13-12-13-2014',...
'BH21-03-05-2015',...
'BR06-03-04-2015',...
'CM17-11-07-2014',...
'ES07-02-23-2015',...
'GR09-07-12-2014',...
'HL17-02-06-2015',...
'J20-08-19-2013',...
'JK15-01-22-2015',...
'JN06-03-28-2015',...
'KT10-10-31-2014',...
'LB24-10-22-2014',...
'LD18-11-07-2014',...
'LW10-06-19-2014',...
'MA16-11-21-2014',...
'N09-12-04-2013',...
'OK09-03-10-2015',...
'OL21-11-23-2014',...
'PC10-02-07-2015',...
'PV22-10-19-2014',...
'R09-09-13-2013',...
'RB23-12-04-2014',...
'SA08-02-22-2015',...
'SM14-03-17-2015',...
'WT06-12-18-2014'}; % ,'RB23-12-04-2014_EyeTracking'  % 'NS24-10-23-2014',... % no kinematics

%% Numbers of behaviors in each behavior list
LastColumnIndex =   [51,...   %A06
                      0,...   %A07
                     83,...   %A18
                      0,...   %AD10
                      0,...   %AG06
                      0,...   %AH10
                    153,...   %AR16
                      0,...   %AS13
                     68,...   %B06
                      0,...   %BG13
                      0,...   %BH21
                      0,...   %BR06
                      0,...   %CM17
                      0,...   %ES07
                    114,...   %GR09 
                      0,...   %HL17
                    144,...   %J20
                      0,...   %JK15
                      0,...   %JN06
                      0,...   %KT10
                      0,...   %LB24
                      0,...   %LD18
                     49,...   %LW10
                      0,...   %MA16
                    104,...   %N09
                      0,...   %OK09
                      0,...   %OL21
                      0,...   %PC10
                      0,...   %PV22
                      0,...   %R09
                    238,...   %RB23
                      0,...   %SA08
                      0,...   %SM14
                      0];     %WT06

%% Selection of infants
InfantID = {'A06','A07','A18','AD10','AG06','AH10','AR16','AS13','B06',...
    'BG13','BH21','BR06','CM17','ES07','GR09','HL17','J20','JK15',...
    'JN06','KT10','LB24','LD18','LW10','MA16','N09','OK09',...
    'OL21','PC10','PV22','R09','RB23','SA08','SM14','WT06'};   %,'RB23_EyeTracking'};   %'NS24', % no kinematics

% infant = menu(['Which infant would you',10,' like to analyze?'],InfantID{:});
% close

for infant = 1:length(InfantID);
%% Initializing and Assigning directory paths....
serverPath1 = ['\\bmi-nas-01\Contreras-UH\Infantdata\Infantdata\Data\',...
    InfantDataAnnotList{infant}];

%% Check if Experiment Start-Stop Indices have been saved already
fullFileName = [serverPath1,'\Behavioral Segmentation\class start-stop times.txt'];
if ~exist(fullFileName, 'file')
  % File does not exist.
  warningMessage = sprintf('Warning: file does not exist:\n%s', fullFileName);
  disp(warningMessage)
  disp('Skipping to next infant data set')
  continue
end

end % repeat for each infant















%% Annotation for Infant Data
InfantDataAnnotList = {'N09-12-04-2013','J20-08-19-2013','B06-10-30-2013',...
    'GR09-07-12-2014','A06-09-28-2013','LW10-06-19-2014','AR16-07-14-2014',...
    'RB23-12-04-2014','A18-08-15-2013'};



InfantID = {'N09','J20','B06','GR09','A06','LW10','AR16','RB23','A18'};
for infant = 1:length(InfantID);
    disp(['Extracting Features using EEG data from ', InfantID{infant}])
%% Initializing and Assigning directory paths....
serverPath1 = ['\\bmi-nas-01\Contreras-UH\Infantdata\Infantdata\code\Zachs_Infant_decoding_files\Data\',InfantDataAnnotList{infant}];
serverPath2 = ['\\172.27.216.40\Contreras-UH\Infantdata\Infantdata\Data\',InfantDataAnnotList{infant}];
desktopPath = ['C:\Users\zrhernan\Infant_decoding_files\Data\',InfantDataAnnotList{infant}];   
cd(desktopPath)

%% Initialization of Data Structures
FEATSELECT = struct('EEG',{},'CLASS',{},'BehaviorSegments',{});

%% Importation of List of Time-Segmented Tasks     
TASKSEGMENTS = importTaskTrialInfo('class start-stop times.txt',2, LastColumnIndex(infant));    % load list of start and stop times for each class
FEATSELECT(1).CLASS.TASKSEGMENTS = TASKSEGMENTS;   % and save into PROCESS structure
disp('____List of tasks imported')

%% Importation of EEG Data
EEG = load('EEGfiles.mat');  % load structure of EEG attributes
disp('____EEG data imported')

%% Resampling Data to 100 Hz
% resampling for EEG data {originally sampled at 1000 Hz}      
    num_EEGchns = size(EEG.uVdata_syncd,2);
    for k=1:num_EEGchns;
        FEATSELECT.EEG.resamp(:,k) = decimate(EEG.uVdata_syncd(:,k),10); 
    end

% resampling for class start-stop times {originally sampled at 1000 Hz}
    DECIMATEDTASKSEGMENTS = FEATSELECT.CLASS.TASKSEGMENTS;
    DECIMATEDTASKSEGMENTS.StartTime = FEATSELECT.CLASS.TASKSEGMENTS.StartTime./10;
    DECIMATEDTASKSEGMENTS.StopTime = FEATSELECT.CLASS.TASKSEGMENTS.StopTime./10;
    FEATSELECT.CLASS.DECIMATEDTASKSEGMENTS = DECIMATEDTASKSEGMENTS; % save into data structure
    
% Initializing time sample sizes
    fs=EEG.srate/10; %[Hz] sampling rate of the decimated signal
    tmax=size(FEATSELECT.EEG.resamp,1)/fs; % number of time samples
    t=0:1/fs:tmax; t1=transpose(t(1:end-1)); % EEG time vector 
    disp(['____ALL data resampled to ' num2str(fs) ' Hz'])

%% Selecting Trials per Behavior
    tasklabels = unique(TASKSEGMENTS.Task); % initialize names given to each behavior
    % Number of seconds lead the onset time (average per class)
    MatchingClassLabels = {'explore', 'imitate',  'observe',  'reach-grasp',  'reach-offer',     'rest'};
    LeadingOnsetperClass =     [5.9,       3.5,        4.5,            1.5,            1.7,        5.4];   % these are in seconds
%     tasklabels = circshift(capitalize(unique(PREPROCESS.CLASS.TASKSEGMENTS.Task)),-1);  % initialize names given to each behavior (alternative
    num_behaviors=max(TASKSEGMENTS.TaskLabel); % number of behaviors to segment
    BehaviorSegments=struct('OnsetTimes',{},'Trials',{}); % create data structure of data segments
    EEGsignalInput = FEATSELECT.EEG.resamp; % initialize EEG signal to use for segmenting by behavior
    % initialize onset times for each behavior
    for cl = 1:num_behaviors; 
        BehaviorSegments(cl).OnsetTimes.Start = DECIMATEDTASKSEGMENTS.StartTime(DECIMATEDTASKSEGMENTS.TaskLabel==cl);
        BehaviorSegments(cl).OnsetTimes.End = DECIMATEDTASKSEGMENTS.StopTime(DECIMATEDTASKSEGMENTS.TaskLabel==cl);
    end  % repeat for all classes per trial

    % Finding trials per class using start and stop onset times
    for cl = 1:num_behaviors;
        disp(['___extracting features from behavior "', MatchingClassLabels{strcmp(MatchingClassLabels,tasklabels(cl))},'"'])
        for p = 1:length(BehaviorSegments(cl).OnsetTimes.Start);
            disp(['_______Trial ', num2str(p)])
            Cstart = BehaviorSegments(cl).OnsetTimes.Start(p);
            Cend = BehaviorSegments(cl).OnsetTimes.End(p);          
        % save all trials separately
            BehaviorSegments(cl).Trials(p).EEGsignal = EEGsignalInput((Cstart:Cend),:);
            oneTrial = BehaviorSegments(cl).Trials(p).EEGsignal;
        % save initialization and completion onsets per trial
            TbeforeOnset = 1*fs; % initialize time before onset to be 1 seconds
            TafterOnset = 3*fs; %LeadingOnsetperClass(strcmp(MatchingClassLabels,tasklabels(cl)))*fs; % initialize time after onset (in seconds)
            if Cstart-TbeforeOnset < 0; continue; end
            BehaviorSegments(cl).Trials(p).OnsetSegment = EEGsignalInput((Cstart-TbeforeOnset:Cstart+TafterOnset),:);
            oneOnsetSegment = BehaviorSegments(cl).Trials(p).OnsetSegment;
            %%% Feature selection
            BehaviorSegments(cl).Trials(p).Features = ComputeFeatures(oneOnsetSegment',fs);
        end     % repeat for all trials for one class
        
        disp('--------------------------------------------------')
        timecompute(toc)
        disp('--------------------------------------------------')
    end  % repeat for all six segmented classes of intended actions
    FEATSELECT.BehaviorSegments = BehaviorSegments;
    disp('____behavioral trials per class selected')
    
    % save to Infantdata path where 'Data' is located
    mkdir(serverPath2,'Feature Selection')
    cd([serverPath2 '/Feature Selection'])
    % create filename
    filename = 'SelectedFeaturesbyTrial.mat';
    save(filename, '-struct', 'FEATSELECT');
    disp('____FEATSELECT data structure saved to lab server folder "Data"')
end