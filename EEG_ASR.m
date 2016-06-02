clc, clear all %, close all
tic  % start computation time

% disp('opening eeglab')
% cd('C:\Users\jgcruz\Documents\MATLAB\eeglab13_4_4b')
% rmpath('C:\Program Files\MATLAB\R2012b\toolbox\eeglab\eeglab12_0_0_0b\ ')
% %eeglab

%% Initialization of Data Structures
PROCESS = struct('EEG_orig',{1},'EEG_new',{1},'EEG_new_trunc',{1});

%% Initializing and Assigning directory paths....
addpath('C:\Users\jgcruz\Documents\MATLAB\eeglab13_4_4b')
addpath('C:\Users\jgcruz\Documents\MATLAB\eeglab13_4_4b\plugins\clean_rawdata0.31')
rmpath('C:\Program Files\MATLAB\R2012b\toolbox\eeglab')

%% Annotation for Infant Data
InfantDataAnnotList = ...
    {'A06-09-28-2013',...
    'A07-08-12-2013',...
    'A18-08-15-2013',...
    'AR16-07-14-2014',...
    'B06-10-30-2013',...
    'BG13-2014-12-13',...
    'CM17-11-07-2014',...
    'GR09-07-12-2014',...
    'J20-08-19-2013',...
    'JK15-01-22-2015',...
    'KT10-10-31-2014',...
    'LB24-10-22-2014',...
    'LD18-11-07-2014',...
    'LW10-06-19-2014',...
    'MA16-11-21-2014',...
    'N09-12-04-2013',...
    'OL21_2014-11-23',...
    'PV22-10-19-2014',...
    'R09-09-13-2013',...
    'RB23-12-04-2014',...
    'WT06-12-18-2014'};

EEGheaderList = ...
    {'autumn-09-28-2013',...
    'A07CC',...
    'asly-08-14-2013',...
    'InfantProject_AR-07-14-2014',...
    'byran-10-30-2013',...
    'InfantProject_000001_BG13',...
    'InfantProject_000001_CM17_2014-11-07',...
    'InfantProject_000001_GR_2014-07-12',...
    'john-08-19-2013',...
    'JK15-01-22-2015'...
    'InfantProject_000001_KT10-10-31-2014',...
    'InfantProject_000001_LB24-2014-10-22-2014',...
    'InfantProject_000001_LD18(2)_2014-11-07',...
    'InfantProject_LW_06-19-14',...
    'K16-11-21-2014','nata-12-04-2013',...
    'InfantProject_000001_OL21_2014-11-23',...
    'InfantProject_PV22-10-19-2014',...
    'rasal-09-13-2013',...
    'InfantProject_RB23',...
    'WT06-12-18-2014'};


%% Selection of infants
% infant = menu(['Which infant would you',10,' like to analyze?'],'N09',...
%     'J20','B06','GR09','A06','LW10','AR16','RB23_EyeTracking','RB23','WT06');
% close
infantID = {'N09','J20','B06','GR09','A06','LW10','AR16'};

%% ASR filter for all infants
tic
for infant = [7:10]
    %% Path directory initializations
    serverPath = ['\\172.27.216.40\Contreras-UH\Infantdata\Infantdata\Data\',InfantDataAnnotList{infant},'\EEG\'];    
    %% Importation of EEG Data
%     PROCESS.EEG_orig = load(['Data\',InfantDataAnnotList{infant},'\EEGfiles.mat']);  % load structure of EEG attributes
%     disp('____EEG data imported')
    %% Opening the VHDR file
    PROCESS.EEG_orig = pop_loadbv(serverPath,[EEGheaderList{infant},'.vhdr']);
    toc
    disp('____EEG data imported')

    OrigEEGdata = PROCESS.EEG_orig;

    %%
    %{
    clean_rawdata(): a wrapper for EEGLAB to call Christian's clean_artifacts.

      Usage:
        >>  EEG = clean_rawdata(EEG, arg_flatline, arg_highpass, arg_channel, arg_noisy, arg_burst, arg_window)

      ------------------ below is from clean_artifacts -----------------------

      This function removes flatline channels, low-frequency drifts, noisy channels, short-time bursts
      and incompletely repaird segments from the data. Tip: Any of the core parameters can also be
      passed in as [] to use the respective default of the underlying functions, or as 'off' to disable
      it entirely.

      Hopefully parameter tuning should be the exception when using this function -- however, there are
      3 parameters governing how aggressively bad channels, bursts, and irrecoverable time windows are
      being removed, plus several detail parameters that only need tuning under special circumstances.

        FlatlineCriterion: Maximum tolerated flatline duration. In seconds. If a channel has a longer
                           flatline than this, it will be considered abnormal. Default: 5

        Highpass :         Transition band for the initial high-pass filter in Hz. This is formatted as
                           [transition-start, transition-end]. Default: [0.25 0.75].

        ChannelCriterion : Minimum channel correlation. If a channel is correlated at less than this
                           value to a reconstruction of it based on other channels, it is considered
                           abnormal in the given time window. This method requires that channel
                           locations are available and roughly correct; otherwise a fallback criterion
                           will be used. (default: 0.85)

        LineNoiseCriterion : If a channel has more line noise relative to its signal than this value, in
                             standard deviations based on the total channel population, it is considered
                             abnormal. (default: 4)

        BurstCriterion : Standard deviation cutoff for removal of bursts (via ASR). Data portions whose
                         variance is larger than this threshold relative to the calibration data are
                         considered missing data and will be removed. The most aggressive value that can
                         be used without losing much EEG is 3. For new users it is recommended to at
                         first visually inspect the difference between the original and cleaned data to
                         get a sense of the removed content at various levels. A quite conservative
                         value is 5. Default: 5.


        WindowCriterion :  Criterion for removing time windows that were not repaired completely. This may
                           happen if the artifact in a window was composed of too many simultaneous
                           uncorrelated sources (for example, extreme movements such as jumps). This is
                           the maximum fraction of contaminated channels that are tolerated in the final
                           output data for each considered window. Generally a lower value makes the
                           criterion more aggressive. Default: 0.25. Reasonable range: 0.05 (very
                           aggressive) to 0.3 (very lax).
    %}

    arg_flatline= -1;
    arg_highpass= [0.05 0.1]; % do you want the HPF to be so high?  I'm not sure what convention is, but I use [0.05 0.1] here...
    arg_channel = -1;
    arg_noisy = -1;
    arg_burst = 3;
    arg_window = -1;
    NewEEGdata = clean_rawdata(OrigEEGdata, arg_flatline, arg_highpass, arg_channel, arg_noisy, arg_burst, arg_window);
    
    toc
    disp(':::::::::::')
    disp('ASR done')
    %% Opening Event Markers (triggers) and plotting them as a bar graph
    eventTriggers = zeros(1,length(EEG.event));
    for trig = 1:length(EEG.event)
        eventTriggers(trig) = EEG.event(trig).latency;
    end
    bar(eventTriggers)
    EEG.eventTriggers = eventTriggers;
    disp('____Event triggers added to data structure')

    %% Input the Experiment Start and Stop Times
    EEG.ExperStartInd = input(['Based on the bar plot, ', 10,...
        'please input the index corresponding to the start of the experiment: ']);
    EEG.ExperEndInd = input(['Based on the bar plot, ', 10,...
        'please input the index corresponding to the end of the experiment: ']);
    disp('____Experiment start and end times added to data structure')

    %% Add modified form of EEG data (where amplitudes are in microVolts)
    EEG.uVdata = double((EEG.data').*(EEG.gain));  % save EEG data into field of 'EEG' structure
    disp('____EEG data in microVolts added to data structure')
    
    %% Sychronizing EEG Data to the Video Recording
    % Shifting Data to match timing of video recording
    % load(['Data\',InfantDataAnnotList{infant},'\ExperStart-EndIndex.txt']); %load start-stop time indices to the whole experiment
    EEG.timesampleStart=EEG.eventTriggers(EEG.ExperStartInd); % sychronized start of EEG data to video recording
    EEG.timesampleEnd=EEG.eventTriggers(EEG.ExperEndInd); % sychronized end of EEG data to video recording
    EEG.uVdata_syncd = EEG.uVdata(EEG.timesampleStart:EEG.timesampleEnd,:); % new EEG set to reflect video recording times
    EEG.uVdata_rmvdCHNS_syncd = EEG.uVdata_rmvdCHNS(EEG.timesampleStart:EEG.timesampleEnd,:); % new EEG set to reflect video recording times
    disp('____video recording-sychronized EEG data added to data structure')

    %% Save clean eeg into data structure
    PROCESS.EEG_ASRv{infant} = NewEEGdata.data;
    PROCESS.EEG_preASRv{infant} = OrigEEGdata.data;
    toc
end
% vis_artifacts(NewEEGdata,OrigEEGdata)
     toc
    disp(':::::::::::')
    disp('visualization done')