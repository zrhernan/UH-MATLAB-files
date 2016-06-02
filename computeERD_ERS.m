
%% Load Trial-Segmented Data from 'Feature Selection' data Structure
load('\\172.27.216.40\Contreras-UH\Infantdata\Infantdata\Data\RB23-12-04-2014\Feature Selection\SelectedFeaturesbyTrial.mat')

%% For calling helper functions
addpath(genpath('\\172.27.216.40\Contreras-UH\Infantdata\Infantdata\code\Zachs_Infant_decoding_files'))

%% Extract onset trials
behavior = 4; % choosing 'reach-to-grasp' behavior
[n_samps, n_chans] = size(BehaviorSegments(behavior).Trials(1).OnsetSegment);
n_trials = length(BehaviorSegments(behavior).Trials);
behavior_onset = zeros(n_samps,n_chans,n_trials);  % preallocate
for trial = 1:n_trials
    behavior_onset(:,:,trial) = BehaviorSegments(behavior).Trials.OnsetSegment;
end

%% Filter onset trials
fs = 100; % sampling frequency
n_f = 3;  % filter order
bpass_freq = [6 9]; % band pass frequencies
behavior_onset_filtered = zeros(n_samps,n_chans,n_trials);  % preallocate
for trial = 1:n_trials
    behavior_onset_filtered(:,:,trial) = filter_data_bpass_NOCELL(behavior_onset(:,:,trial), fs, n_f, bpass_freq);
end

%% Square onset trials
behavior_onset_squared = zeros(n_samps,n_chans,n_trials);  % preallocate
for trial = 1:n_trials
    for chan = 1:n_chans
        behavior_onset_squared(:,chan,trial) = behavior_onset_filtered(:,chan,trial).^2;
    end
end

%% Average onset trials
behavior_onset_mean = mean( behavior_onset_squared, 3 );
behavior_onset_stdev = std( behavior_onset_squared, 0, 3 );

%% Calculate Inter-trial Variance
behavior_onset_var = var( behavior_onset_filtered, 0, 3 );

%% Calculate ERD/ERS
R_bounds = (1:100);

A = behavior_onset_var;

R_avg = mean( A(R_bounds,:),1);
R = repmat(R_avg,n_samps,1);

ERD = ( ( A - R ) ./ R ) * 100;

%% Remove bad channels
cd('C:\Users\zrhernan\Infant_decoding_files\')
load('BrainVision_1020_64ChanLocs.mat'); %files of electrode locations (in Cartesian coordinates)
load('peripheralchnlist.mat');  % list of electrodes along the periphery of the head
load('Data\RB23-12-04-2014\highZchnlist.mat');  % list of electrodes with very high impedances
removechnlist = unique({highZchnlist{:}, peripheralchnlist{:}});  % combine these two lists together

% Save a cell list of channel names
channelOrder = cell(1,length(chanLocs));
for chn=1:length(chanLocs); 
    channelOrder(chn) = cellstr(chanLocs(chn).labels); 
end % save a cell list of channel names
[ channelOrder_rmvdCHNS, ERD_rmvdCHNS ] = removeEEGchannels( channelOrder, removechnlist, 'EEG data', ERD );

%% Extract Central Channels
CentralChans = {'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'Cz',...
                'CP1','CP2','CP3','CP4','CP5','CP6','CPz'};
CentralChansERD = zeros(n_samps,length(CentralChans));            
for Cchan = 1:length(CentralChans)
    if any(strcmp( channelOrder_rmvdCHNS,CentralChans{Cchan} ) ~= 0)
        CentralChansERD(:,Cchan) = ERD_rmvdCHNS(:, strcmp( channelOrder_rmvdCHNS,CentralChans{Cchan} ));
    else
        CentralChansERD(:,Cchan) = nan(n_samps,1);
        continue
    end
end

%% Plot Central Electrodes
figure(1);
h1 = plot((CentralChansERD),'LineWidth',2);
set(gca,'XTick',(1:100:401),'XTickLabel',(-1:3),'xlim',[1 401],'ylim',[500 -500],'YTickLabel',(-500:100:500));
legend(h1, CentralChans{:})

%% Plot Central Electrodes as heat map
figure(2);
imagesc(CentralChansERD');
caxis([-500 500]);      colorbar;
set(gca,'XTick',(1:100:401),'XTickLabel',(-1:3),'xlim',[1 401],...
    'YTick',(1:length(CentralChans)),'YTickLabel',CentralChans,...
    'FontSize',14,'FontWeight','bold');
xlabel('Time (seconds)','FontSize',16,'FontWeight','bold')
ylabel('Central and CentroParietal Channels','FontSize',16,'FontWeight','bold')
title('Percent ERD/ERS for ''Reach-to-Grasp'' Behavior for a 23-Month Boy','FontSize',22,'FontWeight','bold')

%% Plot All Electrodes as heat map
figure(3);
imagesc(ERD_rmvdCHNS');
caxis([-500 500]);      colorbar;
set(gca,'XTick',(1:100:401),'XTickLabel',(-1:3),'xlim',[1 401],...
    'YTick',(1:length(channelOrder_rmvdCHNS)),'YTickLabel',channelOrder_rmvdCHNS,...
    'FontSize',14,'FontWeight','bold');
xlabel('Time (seconds)','FontSize',16,'FontWeight','bold')
ylabel('EEG Channels','FontSize',16,'FontWeight','bold')
title('Percent ERD/ERS for ''Reach-to-Grasp'' Behavior for a 23-Month Boy','FontSize',22,'FontWeight','bold')