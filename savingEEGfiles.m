%% Extraction of EEG data from the Brain Vision VHDR files
%   Other data added to the data structure:
%   >> 'gain' = 0.1 (needed to convert resolution bits to microVolts)
%   >> 'eventTriggers' = (all triggers pressed by the pushbutton during
%       experiment session)
%   >> 'ExperStartInd' = (event marker that corresponds to the initial
%       start of the experiment)
%   >> 'ExperEndInd' = (event marker that corresponds to the initial
%       end of the experiment)
%   >> 'ASRdata' = (EEG data after using ASR method for artifact removal)
%   >> 'uVdata' = (EEG data in microVolts)
%   >> 'channelOrder' = (ordered list of EEG channels corresponding to EEG
%       data)
%   >> 'uVdata_rmvdCHNS' = (EEG data with 'bad channels' removed)
%   >> 'channelOrder_rmvdCHNS' = (list of EEG channels with 'bad channels'
%       removed)
%   >> 'timesampleStart' = (time sample indicating the start of the
%       experiment)
%   >> 'timesampleEnd' = (time sample indicating the end of the experiment)
%   >> 'data_syncd' = (EEG data truncated to the experiment start-end time
%       samples)
clc, clear all, close all

% for calling helper functions
addpath(genpath('\\bmi-nas-01\Contreras-UH\Infantdata\Infantdata\code\Zachs_Infant_decoding_files')) 
%% Select the Extraction of EEG data from the Brain Vision VHDR files
% extractdecision = menu(['Would you like to extract EEG data ', 10,...
%     'from the .vhdr Brain Vision files?'], 'YES', 'NO');
% close
tic;

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
        disp(['Infant folder recognized: ',InfantDataAnnotList{cnt}])
        cnt = cnt + 1;
    end
end
% 
% InfantDataAnnotList = {
% 'A06-09-28-2013',...
% 'A07-08-12-2013',...
% 'A18-08-15-2013',...
% 'AD10-03-24-2015',...
% 'AG06-03-02-2015',...
% 'AH10-04-22-2015',...
% 'AR16-07-14-2014',...
% 'AS13-01-31-2015',...
% 'B06-10-30-2013',...
% 'BG13-12-13-2014',...
% 'BH21-03-05-2015',...
% 'BR06-03-04-2015',...
% 'CM17-11-07-2014',...
% 'ES07-02-23-2015',...
% 'GR09-07-12-2014',...
% 'HL17-02-06-2015',...
% 'J20-08-19-2013',...
% 'JK15-01-22-2015',...
% 'JN06-03-28-2015',...
% 'KT10-10-31-2014',...
% 'LB24-10-22-2014',...
% 'LD18-11-07-2014',...
% 'LW10-06-19-2014',...
% 'MA16-11-21-2014',...
% 'N09-12-04-2013',...
% 'OK09-03-10-2015',...
% 'OL21-11-23-2014',...
% 'PC10-02-07-2015',...
% 'PV22-10-19-2014',...
% 'R09-09-13-2013',...
% 'RB23-12-04-2014',...
% 'SA08-02-22-2015',...
% 'SM14-03-17-2015',...
% 'WT06-12-18-2014'}; % ,'RB23-12-04-2014_EyeTracking'  % 'NS24-10-23-2014',... % no kinematics

%% Selection of infants
% InfantID = {'A06','A07','A18','AD10','AG06','AH10','AR16','AS13','B06',...
%     'BG13','BH21','BR06','CM17','ES07','GR09','HL17','J20','JK15',...
%     'JN06','KT10','LB24','LD18','LW10','MA16','N09','OK09',...
%     'OL21','PC10','PV22','R09','RB23','SA08','SM14','WT06'};   %,'RB23_EyeTracking'};   %'NS24', % no kinematics

% infant = menu(['Which infant would you',10,' like to analyze?'],InfantID{:});
% close

% ExperStartIndList = [3,5,3,7,10,5,7,6,6,11,4];
% ExperEndIndList = [50,43,35,25,29,21,43,16,18,39,24];


for infant = 12:length(InfantDataAnnotList);
    disp(['Acquiring EEG data from ', InfantID{infant}])
%% Path directory initializations
serverPath1 = ['\\172.27.216.40\Contreras-UH\Infantdata\Infantdata\Data\',...
    InfantDataAnnotList{infant},'\EEG'];
serverPath2 = ['\\172.27.216.40\Contreras-UH\Infantdata\Infantdata\code\Zachs_Infant_decoding_files\Data\',...
    InfantDataAnnotList{infant}];
desktopPath = ['C:\Users\zrhernan\Infant_decoding_files\Data\',...
    InfantDataAnnotList{infant}];

%% Check if EEG '.mat' files have been saved already
%
fullFileName = [serverPath1,'\EEGfiles.mat'];
if exist(fullFileName, 'file')
  % File does not exist.
  warningMessage = sprintf('Warning: "EEGfiles.mat" exists:\n%s', fullFileName);
  disp(warningMessage)
  disp('Skipping to next infant data set')
  continue
end
%}
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
%}
%% Opening the VHDR file
EEG = pop_loadbv(serverPath1,[InfantDataAnnotList{infant},'.vhdr']);
EEG.gain = 0.1;      % add gain for changing amplitude order of magnitude
disp('____Brain Vision Recorder EEG file added to data structure')
disp('--------------------------------------------------')
timecompute(toc)
disp('--------------------------------------------------')

%% Using Artifact Subspace Reconstruction (ASR) to remove artifacts
% adding a file path for the 'clean_rawdata' folder
addpath('C:\Program Files\MATLAB\R2012b\toolbox\eeglab\plugins\clean_rawdata0_31')

%     FlatlineCriterion : Maximum tolerated flatline duration. In seconds. If a channel has a longer
arg_flatline= -1;      %  flatline than this, it will be considered abnormal. Default: 5 
        
%             Highpass : Transition band for the initial high-pass filter in Hz. This is formatted as
arg_highpass= [0.05 0.1]; % [transition-start, transition-end]. Default: [0.25 0.75].
         
%     ChannelCriterion : Minimum channel correlation. If a channel is correlated at less than this
arg_channel = -1;     %  value to a reconstruction of it based on other channels, it is considered
%                        abnormal in the given time window. This method requires that channel
%                        locations are available and roughly correct; otherwise a fallback criterion
%                        will be used. (default: 0.85)
% 
%     LineNoiseCriterion : If a channel has more line noise relative to its signal than this value, in
arg_noisy = -1;         %  standard deviations based on the total channel population, it is considered
%                          abnormal. (default: 4)
% 
%     BurstCriterion   :   Standard deviation cutoff for removal of bursts (via ASR). Data portions whose
arg_burst = 3;          %  variance is larger than this threshold relative to the calibration data are
%                          considered missing data and will be removed. The most aggressive value that can
%                          be used without losing much EEG is 3. For new users it is recommended to at
%                          first visually inspect the difference between the original and cleaned data to
%                          get a sense of the removed content at various levels. A quite conservative
%                          value is 5. Default: 5.
% 
% 
%     WindowCriterion :  Criterion for removing time windows that were not repaired completely. This may
arg_window = -1;      %  happen if the artifact in a window was composed of too many simultaneous
%                        uncorrelated sources (for example, extreme movements such as jumps). This is
%                        the maximum fraction of contaminated channels that are tolerated in the final
%                        output data for each considered window. Generally a lower value makes the
%                        criterion more aggressive. Default: 0.25. Reasonable range: 0.05 (very
%                        aggressive) to 0.3 (very lax).  
EEG_temp = clean_rawdata(EEG, arg_flatline, arg_highpass, arg_channel, arg_noisy, arg_burst, arg_window);   
EEG.ASRdata = EEG_temp.data';
disp('____EEG data in microVolts added to data structure')
disp('--------------------------------------------------')
timecompute(toc)
disp('--------------------------------------------------')

%% Import the Experiment Start and Stop Times
load([serverPath1,'\ExperStartStopIndices.mat'])
load([serverPath1,'\EEGTriggerTimes.mat'])
if ExperStartInd == 0 || ExperEndInd == 0
    disp(['Infant ',InfantID{infant},' contains no trigger information'])
    % Saving the resulting variables as a MAT file
    fields_to_save = {'comments', 'nbchan', 'trials',...
    'pnts', 'srate', 'xmin', 'xmax', 'times', 'data', 'chanlocs',...
    'chaninfo', 'ref', 'event', 'eventdescription', 'reject', 'gain',...
    'eventTriggers', 'ExperStartInd', 'ExperEndInd', 'ASRdata', 'uVdata',...
    'channelOrder', 'uVdata_rmvdCHNS', 'channelOrder_rmvdCHNS',...
    'timesampleStart', 'timesampleEnd'};

    cd(serverPath1) % save to Infantdata path where 'Data' is located
    save('EEGfiles.mat', '-struct', 'EEG', fields_to_save{:});
    disp('____EEG data structure saved to lab server folder "Data"')
    clearvars -except serverPath1 InfantDataAnnotList InfantID  infant
    close all
    continue
end

%% Add modified form of EEG data (where amplitudes are in microVolts)
EEG.uVdata = double((EEG.ASRdata).*(EEG.gain));  % save EEG data into field of 'EEG' structure
disp('____EEG data in microVolts added to data structure')

%% Importation of EEG Channels
load('BrainVision_1020_64ChanLocs.mat'); %files of electrode locations (in Cartesian coordinates)
load('peripheralchnlist.mat');  % list of electrodes along the periphery of the head

impedancefileStart = [serverPath1,'\',InfantDataAnnotList{infant},'-start.txt'];
impedancefileEnd = [serverPath1,'\',InfantDataAnnotList{infant},'-end.txt'];
impedancefilePaths = {impedancefileStart; impedancefileEnd};
highZchnlist = createhighZchnlist(impedancefilePaths); % list of electrodes with very high impedances 

removechnlist = transpose(unique([highZchnlist, peripheralchnlist]));  % combine these two lists together

% Save a cell list of channel names
for chn=1:length(chanLocs); 
    EEG.channelOrder(chn) = cellstr(chanLocs(chn).labels); 
end % save a cell list of channel names

[ EEG.channelOrder_rmvdCHNS, EEG.uVdata_rmvdCHNS ] = removeEEGchannels( EEG.channelOrder, removechnlist, 'EEG data', EEG.uVdata );

%{
% Use the 'removechnlist' to remove channels from channel order
chnsremoved=zeros(length(removechnlist),1);
for rmv = 1:length(removechnlist); 
    chnsremoved(rmv)=find(strcmp(EEG.channelOrder,removechnlist{rmv})==1); 
end

EEG.uVdata_rmvdCHNS = EEG.uVdata;        % initialize EEG dataset of all channels
EEG.uVdata_rmvdCHNS(:,chnsremoved)=[];   % and remove undesirable channels from the EEG dataset
EEG.channelOrder_rmvdCHNS = EEG.channelOrder;  % initialize channel order list of all channels
EEG.channelOrder_rmvdCHNS(:,chnsremoved)=[];   % and remove undesirable channels from the EEG channels list
%}
disp('____EEG channels (with "bad" channels removed) added to data structure')
disp('--------------------------------------------------')
timecompute(toc)
disp('--------------------------------------------------')

%% Sychronizing EEG Data to the Video Recording
% Shifting Data to match timing of video recording
EEG.timesampleStart=eventTriggers(ExperStartInd); % sychronized start of EEG data to video recording
EEG.timesampleEnd=eventTriggers(ExperEndInd); % sychronized end of EEG data to video recording
EEG.uVdata_syncd = EEG.uVdata(EEG.timesampleStart:EEG.timesampleEnd,:); % new EEG set to reflect video recording times
EEG.uVdata_rmvdCHNS_syncd = EEG.uVdata_rmvdCHNS(EEG.timesampleStart:EEG.timesampleEnd,:); % new EEG set to reflect video recording times
disp('____video recording-sychronized EEG data added to data structure')
disp('--------------------------------------------------')
timecompute(toc)
disp('--------------------------------------------------')

%% Saving the resulting variables as a MAT file
fields_to_save = {'comments', 'nbchan', 'trials',...
    'pnts', 'srate', 'xmin', 'xmax', 'times', 'data', 'chanlocs',...
    'chaninfo', 'ref', 'event', 'eventdescription', 'reject', 'gain',...
    'ASRdata', 'uVdata','channelOrder', 'uVdata_rmvdCHNS',...
     'channelOrder_rmvdCHNS','timesampleStart', 'timesampleEnd',...
     'uVdata_syncd', 'uVdata_rmvdCHNS_syncd'};

cd(serverPath1) % save to Infantdata path where 'Data' is located
save('EEGfiles.mat', '-struct', 'EEG', fields_to_save{:});
disp('____EEG data structure saved to lab server folder "Data"')
disp('--------------------------------------------------')
timecompute(toc)
disp('--------------------------------------------------')
% cd(serverPath2) % save to Infantdata path where 'Zachs_Infant_decoding_files' are located
% save('EEGfiles.mat', '-struct', 'EEG', fields_to_save{:});
% disp('____EEG data structure saved to lab server folder "Zachs_Infant_decoding_files"')
% disp('--------------------------------------------------')
% timecompute(toc)
% disp('--------------------------------------------------')
% cd(desktopPath) % save to path on Zach's desktop
% save('EEGfiles.mat', '-struct', 'EEG', fields_to_save{:});
% disp('____EEG data structure saved to local Desktop')
% disp('--------------------------------------------------')
% timecompute(toc)
% disp('--------------------------------------------------')
end