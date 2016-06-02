clc, clear all, format compact, %close all
tic  % start computation time

%% Initializing and Assigning directory paths....
serverPath = '\\172.27.216.40\Contreras-UH\Infantdata\Infantdata\Data\N09-12-04-2013';
desktopPath = 'C:\Users\zrhernan\Infant_decoding_files\Data\N09-12-04-2013';
addpath('\\172.27.216.40\Contreras-UH\Infantdata\Infantdata\code\Zachs_Infant_decoding_files') % for calling helper functions
cd(desktopPath)

%% Importation of EEG Data
load([desktopPath, '\EEGfiles.mat']);  % load structure of EEG attributes
disp('____EEG data imported')
disp('--------------------------------------------------')
timecompute(toc)
disp('--------------------------------------------------')

%% Importation of Class Labeling Vector
load([desktopPath, '\classlabel_Observe&Imitate.mat']);     % load class vector
disp('____Target vector imported')

%% Importation of List of Time-Segmented Tasks
TASKSEGMENTS = importTaskTrialInfo([desktopPath,'\Observe&Imitate_start-stop times.txt'],...
    2, 39);  % load list of start and stop times for each 'observe' and 'imitate' class
disp('____List of tasks imported')

%% More Initializations
% initialize number of behaviors to classify
num_behaviors = max(unique(classlabel));

% initialize channel and time sample number of EEG data
[num_samples, num_EEGchns] = size(uVdata_syncd);

% store names for each behavior
BehaviorName = flipud(unique(TASKSEGMENTS.Task)); 

%% Resampling Data to 100 Hz
% resampling for EEG data {originally sampled at 1000 Hz}      
fs = srate/10; %[Hz] sampling rate of the decimated signal
EEG_resamp = zeros(ceil(num_samples/10), num_EEGchns);
for k = 1:num_EEGchns;
    EEG_resamp(:,k) = decimate(uVdata_syncd(:,k),10); 
end
disp(['____ALL data resampled to ' num2str(fs) ' Hz'])

%% Resample class start-stop times {originally sampled at 1000 Hz}
NEWTASKSEGMENTS = TASKSEGMENTS;
NEWTASKSEGMENTS.StartTime = TASKSEGMENTS.StartTime./10;
NEWTASKSEGMENTS.StopTime = TASKSEGMENTS.StopTime./10;
    
%% Selecting Trials per Behavior
Behaviors=struct('OnsetTimes',{},'Trials',{});

% initialize onset times for each behavior
for cl = 1:num_behaviors; 
    Behaviors(cl).OnsetTimes.Start = NEWTASKSEGMENTS.StartTime(NEWTASKSEGMENTS.TaskLabel==cl);
    Behaviors(cl).OnsetTimes.End = NEWTASKSEGMENTS.StopTime(NEWTASKSEGMENTS.TaskLabel==cl);
end  % repeat for all behaviors per trial

chns = (1:num_EEGchns);
chnCombos = nchoosek(chns,2);       scales = (1:512); wname = 'morl';
chnCombos1 = chnCombos(:,1);        chnCombos2 = chnCombos(:,2);
% Finding trials per behavior using start and stop onset times
for cl = 1:num_behaviors;
    for p = 1:length(Behaviors(cl).OnsetTimes.Start);
        Cstart = Behaviors(cl).OnsetTimes.Start(p);
        Cend = Behaviors(cl).OnsetTimes.End(p);

    % save all trials separately
        Behaviors(cl).Trials(p).EEGsignal = EEG_resamp((Cstart:Cend),:);
        oneTrial = Behaviors(cl).Trials(p).EEGsignal;
    %    
    % compute wavelet coherence
        EEGCombins1 = oneTrial(:,chnCombos1);
        EEGCombins2 = oneTrial(:,chnCombos2);
        %
%             if matlabpool('size') == 0; matlabpool('open'); end
        for combin = 1:length(chnCombos);
            disp('Performing wavelet coherence')
            disp([num2str(sRound((combin/length(chnCombos))*100,2)),...
                ' Percent of EEG Combinations finished'])                
            disp(['...between ',channelOrder{chns(chnCombos(combin,1))}, ' & ',...
            channelOrder{chns(chnCombos(combin,2))}, ' channels'])
            disp(['For ',BehaviorName{cl},' & Trial ',num2str(p)])

            % initialize the EEG channels pairs to use for classification
            eeg1 = EEGCombins1(:,combin);
            eeg2 = EEGCombins2(:,combin);

            % apply coherence and cross spectrum to each trial per class
            [wcoh(:,:,combin), wcs(:,:,combin)] = wcoher(eeg1,eeg2,scales,wname);
            wcs_modulus = abs(wcs(:,:,combin));
            [maxwcs(combin), max_ind] = max(wcs_modulus(:));
            [maxscal(combin),maxtimesamp(combin)] = find(wcs_modulus == wcs_modulus(max_ind));
        end % repeat for all paired channel combinations
                
        % save variables into data structure 'Behaviors'
        Behaviors(cl).Trials(p).maxwcs = maxwcs;
        Behaviors(cl).Trials(p).wcoh = wcoh;
        Behaviors(cl).Trials(p).wcs = wcs;
        Behaviors(cl).Trials(p).maxscal = maxscal;
        Behaviors(cl).Trials(p).maxtimesamp = maxtimesamp;
        
        % then clear out variables for next round of EEG channel pair combinations
        clear wcoh wcs wcs_modulus maxwcs maxscal maxtimesamp
        disp('--------------------------------------------------')
        timecompute(toc)
        disp('--------------------------------------------------')
        %}
    end     % repeat for all trials for one behavior
    disp('--------------------------------------------------')
    timecompute(toc)
    disp('--------------------------------------------------')
end  % repeat for all behaviors
disp('____all trials per behavior stored and features computed')
disp('--------------------------------------------------')
timecompute(toc)
disp('--------------------------------------------------')
% matlabpool('close')
return
%% Extract trial features
% Extract an array of scalar features by trial and concatenate to form a 
% [trial-by-channel] matrix with concatenated classes

Target = zeros(1, length(NEWTASKSEGMENTS)); % preallocate
MaxWSCPerTrial = zeros(length(NEWTASKSEGMENTS), num_EEGchns); % preallocate
tnum = 0;

% extract the trials from the 'Behaviors' data structure
for cl = 1:num_behaviors;
    for p = 1:length(Behaviors(cl).Trials);
        tnum = tnum + 1;
        MaxWSCPerTrial(tnum,:) = Behaviors(cl).Trials(p).maxwcs;
        Target(tnum) = cl;
    end
end


% Plot the features as a scatter plot (for demonstration)
%
fig_hndl=figure;
gscatter(MaxWSCPerTrial(:,1),MaxWSCPerTrial(:,6),Target,'rb','v^')

% additional parameters to enhance the plot
xlabel(channelOrder{chns(chnCombos(combin,1))},'FontSize',12,'FontWeight','bold')
ylabel(channelOrder{chns(chnCombos(combin,2))},'FontSize',12,'FontWeight','bold')
set(gca,'FontWeight','bold','FontSize',12)

enhanceFigure(fig_hndl); % enhance figure for export
mkdir('Plots') % make new 'Plots' folder
filename = 'scatterplotofmaxWSC'; % name of saved file
saveas(gcf, [filename, '.fig'], 'fig'); % save as .fig file
hgexport(gcf, [filename, '.tif'], hgexport('factorystyle'), 'Format', 'tiff'); % save as .tif file       
cd('..') % go back to main folder
%}
disp('--------------------------------------------------')
timecompute(toc)
disp('--------------------------------------------------')