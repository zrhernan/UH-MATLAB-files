clc, clear all, format compact, %close all
tic  % start computation time

%% for calling helper functions
addpath(genpath('\\172.27.216.40\Contreras-UH\Infantdata\Infantdata\code\Zachs_Infant_decoding_files'))

%% Opening directory path for Prasad's position files
cd('yourdirectorypath')

%% Import Kinematics Data
% for importing position files
    FullSetKine = load('yourfilehere.mat');

%% Import List of Time-Segmented Tasks
    TASKSEGMENTS = importTaskTrialInfo('class start-stop times.txt');
    disp('____List of tasks imported')    

%% Resample Data to 100 Hz
    fs = 100; %[Hz] sampling rate of the decimated signal
    
    % resample for gravity-compensated magnitude acceleration (GCMA) data {originally sampled at 128 Hz}
    num_KINEsnsrs = size(FullSetKine,2);
    for p = 1:num_KINEsnsrs;
        FullSetAccel_resamp = resample(FullSetKine,100,128); 
    end

    % decimate class start-stop times {originally sampled at 1000 Hz}
    DECIMATEDTASKSEGMENTS = TASKSEGMENTS;
    DECIMATEDTASKSEGMENTS.StartTime = ceil(TASKSEGMENTS.StartTime./10);
    DECIMATEDTASKSEGMENTS.StopTime = ceil(TASKSEGMENTS.StopTime./10);
    disp(['____ALL data resampled to ' num2str(fs) ' Hz'])

%% Select Trials per Behavior
    tasklabels = unique(TASKSEGMENTS.Task); % initialize names given to each behavior
    % Number of seconds lead the onset time (average per class)
    MatchingClassLabels = {'explore', 'imitate',  'observe',  'reach-grasp',  'reach-offer',     'rest'};
    LeadingOnsetperClass =     [5.9,       3.5,        4.5,            1.5,            1.7,        5.4];   % these are in seconds
    num_behaviors=max(TASKSEGMENTS.TaskLabel); % number of behaviors to segment
    BehaviorSegments=struct('OnsetTimes',{},'Trials',{}); % create data structure of data segments
    signalInput = FullSetAccel_resamp; % initialize signal to use for segmenting by behavior
    
    % initialize onset times for each behavior
    for cl = 1:num_behaviors;
        classidx = DECIMATEDTASKSEGMENTS.TaskLabel==cl;
        BehaviorSegments(cl).OnsetTimes.Start = DECIMATEDTASKSEGMENTS.StartTime(classidx);
        BehaviorSegments(cl).OnsetTimes.End = DECIMATEDTASKSEGMENTS.StopTime(classidx);
    end  % repeat for all classes per trial

    % Finding trials per class using start and stop onset times
    for cl = 1:num_behaviors;
        BehaviorName = MatchingClassLabels{strcmp(MatchingClassLabels,tasklabels(cl))};
        BehaviorSegments(cl).Name = BehaviorName; %add to data structure
        disp(['___extracting features from behavior "', BehaviorName,'"'])
        for p = 1:length(BehaviorSegments(cl).OnsetTimes.Start);
            disp(['_______Trial ', num2str(p)])
            BehaviorSegments(cl).Trials(p).TrialNumber = p;
            Cstart = BehaviorSegments(cl).OnsetTimes.Start(p);
            Cend = BehaviorSegments(cl).OnsetTimes.End(p);
        % save all trials separately
            BehaviorSegments(cl).Trials(p).Accelsignal = signalInput((Cstart:Cend),:);
            oneTrial = BehaviorSegments(cl).Trials(p).Accelsignal;
        % save initialization and completion onsets per trial
            TbeforeOnset = 1*fs; % initialize time before onset to be 1 seconds
            TafterOnset = 3*fs; %LeadingOnsetperClass(strcmp(MatchingClassLabels,tasklabels(cl)))*fs; % initialize time after onset (in seconds)
            if Cstart-TbeforeOnset < 0; continue; end
            BehaviorSegments(cl).Trials(p).OnsetSegment = signalInput((Cstart-TbeforeOnset:Cstart+TafterOnset),:);
            oneOnsetSegment = BehaviorSegments(cl).Trials(p).OnsetSegment;
%=========================================================================%
        % Feature selection
            BehaviorSegments(cl).Trials(p).Features = ComputeFeatures(oneOnsetSegment',fs);
%=========================================================================%
        % compute histogram
            histx=(-3:0.01:3);   % initialize
            BehaviorSegments(cl).Trials(p).histogram = hist(oneTrial,histx);
%=========================================================================%          
        end     % repeat for all trials for one class
        disp('--------------------------------------------------')
        timecompute(toc)
        disp('--------------------------------------------------') 
    end  % repeat for all six segmented classes of intended actions
    disp('____behavioral trials per class selected')