clc, clear all, format compact, %close all
tic  % start computation time

%% for calling helper functions
addpath(genpath('\\172.27.216.40\Contreras-UH\Infantdata\Infantdata\code\Zachs_Infant_decoding_files'))

%% remove this path (for 'resample' function)
rmpath(genpath('C:\Program Files\MATLAB\R2012b\toolbox\signal\signal\ja'))

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

%% Selection of infants
% infant = menu(['Which infant would you',10,' like to analyze?'],InfantID{:});
% close

%% List of IMU sensor types
sensorIDList_post2013 = {'SI-000719','SI-000775','SI-000722','SI-000773',...
    'SI-000738','SI-000708'}; % OPAL sensor IDs for experiments after 2013
sensorIDList_2013 = {'SI-000708','SI-000719','SI-000722','SI-000738',...
    'SI-000742','SI-000773'}; % OPAL sensor IDs for 2013 experiments
sensorBPList = { 'ExperimentersLeftArm','ExperimentersRightArm',...
    'InfantsLeftArm','InfantsRightArm','InfantsTrunk','InfantsForehead'};

for infant = 5%1:length(InfantID);
    if infant == 2; continue; end
    if ~strcmp(InfantID,'N09'); continue; end
%% Initializing and Assigning directory paths....
    serverPath1 = ['\\bmi-nas-01\Contreras-UH\Infantdata\Infantdata\Data\',...
        InfantDataAnnotList{infant}];

%% Check if Experiment Start-Stop Times have been saved already
    fullFileName = [serverPath1,'\Behavioral Segmentation\class start-stop times recheck.txt'];
    if ~exist(fullFileName, 'file')
      % File does not exist.
      warningMessage = sprintf('Warning: file does not exist:\n%s', fullFileName);
      disp(warningMessage)
      disp('Skipping to next infant data set')
      continue
    end
    
%% Check if Kinematics Features have been saved already
    fullFileName = [serverPath1,'\Feature Selection\SelectedAccelFeaturesbyTrial.mat'];
    if exist(fullFileName, 'file')
      % File does not exist.
      warningMessage = sprintf('Warning: Kinematics Features file already exists:\n%s', fullFileName);
      disp(warningMessage)
      disp('Skipping to next infant data set')
      continue
    end    
    disp(['Extracting Features using EEG or kinematics data from ', InfantID{infant}])

%% Import List of Time-Segmented Tasks
    TASKSEGMENTS = importTaskTrialInfo([serverPath1,'\Behavioral Segmentation\class start-stop times recheck.txt']);%, LastColumnIndex(infant));    % load list of start and stop times for each class
    disp('____List of tasks imported')

%% Specify which sensor ID set to use (depends on the year experiment was conducted)
    if ~isempty(strfind(InfantDataAnnotList{infant},'2014')) || ~isempty(strfind(InfantDataAnnotList{infant},'2015'));
        sensorIDList = sensorIDList_post2013;
    elseif ~isempty(strfind(InfantDataAnnotList{infant},'2013'));
        sensorIDList = sensorIDList_2013;
    end

%% Import Kinematics Data
    extractedKINEFiles = {'Acc_n_syncd','GCMA_syncd','vel_filt_syncd'};
    KINE = struct(extractedKINEFiles{1},{},extractedKINEFiles{2},{},extractedKINEFiles{3},{});
    for sens = 1:6
        fullFileName_KINE = [serverPath1,'\Kinematics\',sensorBPList{sens},'_',sensorIDList{sens},'.mat'];
        if ~exist(fullFileName_KINE, 'file')
        % IMU sensor data does not exist.
            warningMessage = sprintf('Warning: IMU sensor data does not exist:\n%s', fullFileName);
            disp(warningMessage)
            disp('Skipping to next IMU sensor')
            continue
        else
            contents = whos('-file',fullFileName_KINE);
            % Check if synchronized kinematics data has been saved already
            if all(ismember(extractedKINEFiles, {contents.name}))
                KINE(sens) = load(fullFileName_KINE, extractedKINEFiles{:});  % load kinematics data
                disp(['____Kinematics data for ',sensorBPList{sens},' Sensor imported'])    
            else
            % Synchronized data does not exist.
                warningMessage = sprintf('Warning: Synchronized acceleration data does not exist:\n%s', fullFileName);
                disp(warningMessage)
                disp('Skipping to next IMU sensor')
                continue
            end
        end
    end
    FullSetKine = [KINE.Acc_n_syncd, KINE.GCMA_syncd, KINE.vel_filt_syncd];

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
    
%%  Band-pass filter within a specific band
%     n_f = 3;  % filter order
%     bpass_freq=[0.001 6]; % band-pass frequencies
%     FullSetAccel_filtered = filter_data_bpass_NOCELL(FullSetAccel_resamp, fs, n_f, bpass_freq);
%     disp('____ALL data filtered')   

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
%             BehaviorSegments(cl).Trials(p).Features = ComputeFeatures(oneOnsetSegment',fs);
%=========================================================================%
        % compute histogram
            histx=(-3:0.01:3);   % initialize
            BehaviorSegments(cl).Trials(p).histogram = hist(oneTrial,histx);
            %{
%=========================================================================%
        % compute minimum value
            BehaviorSegments(cl).Trials(p).min = min(oneTrial);
%=========================================================================%
        % compute maximum value
            BehaviorSegments(cl).Trials(p).min = max(oneTrial);
%=========================================================================%
        % compute standard deviation value
            BehaviorSegments(cl).Trials(p).std = std(oneTrial);
%=========================================================================%
        % compute kurtosis value
            BehaviorSegments(cl).Trials(p).kurt = kurtosis(oneTrial);
%=========================================================================%
        % compute Shannon's entropy value
            BehaviorSegments(cl).Trials(p).kurt = wentropy(oneTrial,'shannon');
%=========================================================================%
        % compute spectral estimation           
            for d=1:num_KINEsnsrs;
            BehaviorSegments(cl).Trials(p).SPECTROGRAM(d).DimensionNumber = d;    
            % Fast Fourier and Thompson's Multitaper Method Initializations
            nw=4; nfft_fft = 512;   nfft_pmtm = 2^nextpow2(size(oneTrial(:,d),1));
            if nfft_pmtm < 100; % keep out extremely small trial sizes
                continue
            end
            
            % calculate the Fourier Transform (using Thompson's Multi-Taper Method)
                [BehaviorSegments(cl).Trials(p).pmtm_psd(:,d),...
                    BehaviorSegments(cl).Trials(p).pmtm_confid(:,:,d),...
                    BehaviorSegments(cl).Trials(p).pmtm_freq(:,d)]...
                    = pmtm(oneTrial(:,d), nw, nfft_pmtm, fs, 0.95);                
                
            % calculate the Fourier Transform (using Fast Fourier Algorithm)    
                BehaviorSegments(cl).Trials(p).fft(:,d) = ...
                    fft(oneTrial(:,d),nfft_fft)/size(oneTrial(:,d),1);

            % calculate the short-time Fast Fourier transform (STFT)
            %{
                wndw =20;     novrlp = wndw-1;      % initialize
                % Initialization Onset transitions    
                    Onset_trial = oneOnsetSegment; init_spect_flg=1;
                % Completion Onset transitions
%                     Onset_trial = oneCompOnset; comp_spect_flg=1;
                OnsetTriallength = length(Onset_trial);  % spectrogram number of frequency samples is the length of the trial period
                nfft_s = 2^nextpow2(OnsetTriallength);
                [ BehaviorSegments(cl).Trials(p).SPECTROGRAM(d).stft,...
                    BehaviorSegments(cl).Trials(p).SPECTROGRAM(d).freq, ...
                    BehaviorSegments(cl).Trials(p).SPECTROGRAM(d).time,...
                    BehaviorSegments(cl).Trials(p).SPECTROGRAM(d).stftpsd ] = ...
                    spectrogram(Onset_trial(:,d), wndw, novrlp, nfft_s, fs); % generate STFT spectrogram
            %}
            end % repeat for all channels
            %}
%=========================================================================%          
        end     % repeat for all trials for one class
%         disp('--------------------------------------------------')
%         timecompute(toc)
%         disp('--------------------------------------------------') 
    end  % repeat for all six segmented classes of intended actions
    disp('____behavioral trials per class selected')

    % save to Infantdata path where 'Data' is located
    mkdir(serverPath1,'Feature Selection')
    cd([serverPath1 '/Feature Selection'])
    % create filename
    filename = 'SelectedAccelFeaturesbyTrial.mat';
    save(filename, 'BehaviorSegments','-v7.3');
    disp('____FEATSELECT data structure saved to lab server folder "Data"')
    close all
    clearvars -except InfantDataAnnotList InfantID sensorIDList_post2013 sensorIDList_2013 sensorBPList 
end  % repeat for each infantID