clc, clear all, format compact, %close all
tic  % start computation time

%% Annotation for Infant Data
InfantDataAnnotList = {'N09-12-04-2013','J20-08-19-2013','B06-10-30-2013'};
infant = menu(['Which infant would you',10,' like to analyze?'],'N09','J20','B06');
close all
%% Initialization of Data Structures
PREPROCESS = struct('EEG',{},'KINE',{},'CLASS',{});
% PREPROCESS(1).KINE = struct('InfantForehead',{},'InfantArmL',{},'InfantArmR',{},'ActorArmL',{},'ActorArmR',{});

%% Initializing and Assigning directory paths....
InfantData_dir = 'Z:\Infantdata\Infantdata\code\Zachs_Infant_decoding_files';
cd(InfantData_dir)

%% Importation of Class Labeling Vector
load(['Data\',InfantDataAnnotList{infant},'\classlabel.mat']);     % load class vector
% load(['Data\',InfantDataAnnotList{infant},'\classlabel_Observe&Imitate.mat']);     % load class vector
PREPROCESS(1).CLASS.classlabel=classlabel; clear classlabel % and save into data structure
disp('____Target vector imported')

%% Importation of List of Time-Segmented Tasks
LastColumnIndex = [104,144,68,114,51,49,153,238,83];              LastColumnIndex_2class = [39,46,19]; %for binary class dataset
LastColumnIndex_7C_J20 = 162; 
TASKSEGMENTS = importTaskTrialInfo(['Data\',InfantDataAnnotList{infant},'\class start-stop times.txt'],...
    2, LastColumnIndex(infant));  %LastColumnIndex_7C_J20); %    % load list of start and stop times for each class
% TASKSEGMENTS = importTaskTrialInfo(['Data\',InfantDataAnnotList{infant},'\Observe&Imitate_start-stop times.txt'],...
%     2, LastColumnIndex_2class);  % load list of start and stop times for each 'observe' and 'imitate' class
PREPROCESS.CLASS.TASKSEGMENTS = TASKSEGMENTS;   % and save into PROCESS structure
disp('____List of tasks imported')

%% Importation of EEG Data
PREPROCESS.EEG = load(['Data\',InfantDataAnnotList{infant},'\EEGfiles.mat']);  % load structure of EEG attributes
% figure; 
% ha=tight_subplot(1,1,[.01 .03],[.1 .01],[.01 .01]);   axes(ha(1))
% imagesc(PROCESS.EEG.uVdata'); title('Raw EEG Data');
% set(gca,'YTick',(1:length(PROCESS.EEG.channelOrder)),...
%     'YTickLabel',PROCESS.EEG.channelOrder,...
%     'XTick',(PROCESS.EEG.times./PROCESS.EEG.srate),...
%     'XTickLabel',(PROCESS.EEG.xmin:100:PROCESS.EEG.xmax))
% colorbar('Position',[0.931 0.11 0.014 0.816])
disp('____EEG data imported')

%% Importation of Kinematics Data
%{
PREPROCESS.KINE.InfantH = load ('Data\',InfantDataAnnotList{infant},'\InfantsForehead_SI-000773.mat'); % save kinematics Data into each body sensor
PREPROCESS.KINE.InfantT = load ('Data\',InfantDataAnnotList{infant},'\InfantsTrunk_SI-000742.mat'); % save kinematics Data into each body sensor
PREPROCESS.KINE.InfantLA = load ('Data\',InfantDataAnnotList{infant},'\InfantsLeftArm_SI-000722.mat'); % save kinematics Data into each body sensor
PREPROCESS.KINE.InfantRA = load ('Data\',InfantDataAnnotList{infant},'\InfantsRightArm_SI-000738.mat'); % save kinematics Data into each body sensor
PREPROCESS.KINE.FullSet = horzcat(PREPROCESS.KINE.InfantH.Acc_n,...
    PREPROCESS.KINE.InfantT.Acc_n,PREPROCESS.KINE.InfantLA.Acc_n,PREPROCESS.KINE.InfantH.Acc_n);
PREPROCESS.KINE.InfantH.GCMA = transpose(PREPROCESS.KINE.InfantH.MA);

% Acquiring Trigger data
 for us=1:length(PREPROCESS.KINE.InfantH.Annotations.Time); 
     PREPROCESS.CLASS.KINETrigEdgeOnset(us,1)=PREPROCESS.KINE.InfantH.Annotations.Time(us)-min(PREPROCESS.KINE.InfantH.Time_stamps);
 end
PREPROCESS.CLASS.KINETrigEdgeOnset = double(PREPROCESS.CLASS.KINETrigEdgeOnset.*(128/1e+6));  %re-format so that onset values are in 128 Hz samples 

% Shifting Data to match timing of video recording
KINEStart=PREPROCESS.CLASS.KINETrigEdgeOnset(7); % sychronized start of kinematics data to video recording
KINEEnd=PREPROCESS.CLASS.KINETrigEdgeOnset(83); % sychronized end of kinematics data to video recording
PREPROCESS.KINE.FullSet_vidsync = PREPROCESS.KINE.FullSet(KINEStart:KINEEnd,:); % new kinematics set to reflect video recording times
PREPROCESS.KINE.InfantH.GCMA_vidsync = PREPROCESS.KINE.InfantH.GCMA(KINEStart:KINEEnd,:); % new kinematics set to reflect video recording times
%}
%% Resampling Data to 100 Hz
% resampling for EEG data {originally sampled at 1000 Hz}      
    num_EEGchns = size(PREPROCESS.EEG.uVdata_syncd,2);
    for k=1:num_EEGchns;
        PREPROCESS.EEG.resamp(:,k) = decimate(PREPROCESS.EEG.uVdata_syncd(:,k),10); 
    end

% resampling for gravity-comensated magnitude acceleration (GCMA) data {originally sampled at 128 Hz}
%     PREPROCESS.KINE.InfantH.resampGCMA = resample(PREPROCESS.KINE.InfantH.GCMA_vidsync,100,128);

% resampling for class start-stop times {originally sampled at 1000 Hz}
    DECIMATEDTASKSEGMENTS = PREPROCESS.CLASS.TASKSEGMENTS;
    DECIMATEDTASKSEGMENTS.StartTime = PREPROCESS.CLASS.TASKSEGMENTS.StartTime./10;
    DECIMATEDTASKSEGMENTS.StopTime = PREPROCESS.CLASS.TASKSEGMENTS.StopTime./10;
    PREPROCESS.CLASS.DECIMATEDTASKSEGMENTS = DECIMATEDTASKSEGMENTS; % save into data structure

% Precise alignment of both data streams
% Take out extra samples in the beginning. Should only be 1-5 samples to remove
%     samplesrmvd = abs(length(PREPROCESS.EEG.resamp) - length(PREPROCESS.KINE.InfantH.resampGCMA));
%     if samplesrmvd > 5;
%         disp('Resampled biomedical data streams larger than five samples, re-check syncing alignment.');
%     else
%         if length(PREPROCESS.EEG.resamp) > length(PREPROCESS.KINE.InfantH.resampGCMA);
%             PREPROCESS.EEG.resamp = PREPROCESS.EEG.resamp(samplesrmvd+1:end,:);
%         elseif length(PREPROCESS.KINE.InfantH.resampGCMA) > length(PREPROCESS.EEG.resamp);
%             PREPROCESS.KINE.InfantH.resampGCMA = PREPROCESS.KINE.InfantH.resampGCMA(samplesrmvd+1:end,:);
%         end
%     end
    
% Initializing time sample sizes
    fs=100; %[Hz] sampling rate of the decimated signal
    tmax=size(PREPROCESS.EEG.resamp,1)/fs; % number of time samples
    t=0:1/fs:tmax; t1=transpose(t(1:end-1)); % EEG time vector
%     t2 = linspace(0,(length(PREPROCESS.KINE.InfantH.resampGCMA(:,1))/100),...
%         length(PREPROCESS.KINE.InfantH.resampGCMA(:,1)));  % Head GCMA time vector  
    disp(['____ALL data resampled to ' num2str(fs) ' Hz'])

%% Transform signal from MAT to CELL
%
% Initialization of cell array
    PREPROCESS.EEG.resamp_cell=cell(1,num_EEGchns);
    
% For creating the cell array of EEG data
    for k=1:num_EEGchns;
        PREPROCESS.EEG.resamp_cell{1,k}=PREPROCESS.EEG.resamp(:,k);
    end
    
%%  Extract Delta Band EEG using a band-pass filter
    n_f=3;  % filter order
    delta_freq=[0.2 32]; % band passing into the delta frequency
    PREPROCESS.EEG.filtered_cell=filter_data_bpass(PREPROCESS.EEG.resamp_cell, fs, n_f, delta_freq);
    PREPROCESS.EEG.filtered = cell2mat(PREPROCESS.EEG.filtered_cell);
     disp('____EEG data filtered')
    
%% Standardization across channels (Z-scores: subtracting mean from each data point and deviding by standard deviation )
%   tip: mean 'should be already' zero after filtering, check this if your filter is working correctly
    PREPROCESS.EEG.zscores_filtered = zscore(PREPROCESS.EEG.filtered);
    disp('____EEG data standardized')

%% Plot sample EEG signal 'CPz'
    %{
    figure(10);     sampCHN = find(strcmp(PREPROCESS.EEG.channelOrder_rmvdCHNS,'CPz'));
% plot raw and filtered EEG on the first subplot
    subplot(2,1,1); plot(t1,PREPROCESS.EEG.resamp_cell{sampCHN},'k',t1,PREPROCESS.EEG.filtered_cell{sampCHN},'r'); % plot raw and filtered CPz
    hold on; BT_handl=plot(t1(1:end-1),(PREPROCESS.CLASS.classlabel*max(PREPROCESS.EEG.filtered_cell{sampCHN}))+min(PREPROCESS.EEG.filtered_cell{sampCHN}),'g--'); % also plot the class boundaries
    %set(get(BT_handl,'BaseLine'),'LineStyle',':');   set(BT_handl,'MarkerFaceColor','red','MarkerSize',0,'LineWidth',2);
    title(['Signal Pre-processing: '  PREPROCESS.EEG.channelOrder_rmvdCHNS{sampCHN} ' electrode'],'FontSize',14,'FontWeight','bold');
    ylabel('Voltage (\muV)','FontSize',12,'FontWeight','bold');
    xlim([min(t1) max(t1)]); ylim([min(PREPROCESS.EEG.raw_vidsync(:,sampCHN)) max(PREPROCESS.EEG.raw_vidsync(:,sampCHN))]);
         legend('raw EEG',['filtered EEG (zero-phase, 3rd order Butterworth, [' delta_freq(1) '-' delta_freq(2) ' Hz])']);
% then plot the standardized, filtered signal in the second subplot
    subplot(2,1,2); plot(t1,PREPROCESS.EEG.zscores_filtered(:,sampCHN),'b');
    xlim([min(t1) max(t1)]); legend('standardized EEG (zero mean, unit variance)'); 
    ylim([min(PREPROCESS.EEG.zscores_filtered(:,sampCHN)) max(PREPROCESS.EEG.zscores_filtered(:,sampCHN))]);
    xlabel('Time (s)','FontSize',12,'FontWeight','bold'); ylabel('Z-Score (\sigma)','FontSize',12,'FontWeight','bold');
    hold on; BT_handl=plot(t1(1:end-1),(PREPROCESS.CLASS.classlabel*max(PREPROCESS.EEG.zscores_filtered(:,sampCHN)))+min(PREPROCESS.EEG.zscores_filtered(:,sampCHN)),'g--'); % also plot the class boundaries
    %}
%% Selecting Trials per Behavior
    tasklabels = flipud(unique(TASKSEGMENTS.Task)); % initialize names given to each behavior
%     tasklabels = circshift(capitalize(unique(PREPROCESS.CLASS.TASKSEGMENTS.Task)),-1);  % initialize names given to each behavior (alternative
    num_behaviors=max(TASKSEGMENTS.TaskLabel); % number of behaviors to segment
    BehaviorSegments=struct('OnsetTimes',{},'Trials',{}); % create data structure of data segments
    EEGsignalInput = PREPROCESS.EEG.resamp; %PREPROCESS.EEG.zscores_filtered; % initialize EEG signal to use for segmenting by behavior
    % initialize onset times for each behavior
    for cl = 1:num_behaviors; 
        BehaviorSegments(cl).OnsetTimes.Start = DECIMATEDTASKSEGMENTS.StartTime(DECIMATEDTASKSEGMENTS.TaskLabel==cl);
        BehaviorSegments(cl).OnsetTimes.End = DECIMATEDTASKSEGMENTS.StopTime(DECIMATEDTASKSEGMENTS.TaskLabel==cl);
    end  % repeat for all classes per trial

    % Finding trials per class using start and stop onset times
    for cl = 1:num_behaviors;
        for p = 1:length(BehaviorSegments(cl).OnsetTimes.Start);
            Cstart = BehaviorSegments(cl).OnsetTimes.Start(p);
            Cend = BehaviorSegments(cl).OnsetTimes.End(p);          
        % save all trials separately
            BehaviorSegments(cl).Trials(p).EEGsignal = EEGsignalInput((Cstart:Cend),:);
            oneTrial = BehaviorSegments(cl).Trials(p).EEGsignal;
        % save initialization and completion onsets per trial
            trans_intrvl = 200; % initialize number of samples before and after trial onset to use for segmenting
            BehaviorSegments(cl).Trials(p).initOnset = EEGsignalInput((Cstart-trans_intrvl:Cstart+trans_intrvl),:);
            BehaviorSegments(cl).Trials(p).compOnset = EEGsignalInput((Cstart-trans_intrvl:Cstart+trans_intrvl),:);
            oneInitOnset = BehaviorSegments(cl).Trials(p).initOnset;
            oneCompOnset = BehaviorSegments(cl).Trials(p).compOnset;
        % compute histogram
            histx=(-3:0.01:3);   % initialize
            BehaviorSegments(cl).Trials(p).histogram = hist(oneTrial,histx);
        % compute minimum value
            BehaviorSegments(cl).Trials(p).min = min(oneTrial);       
        % compute spectral estimation           
            for chn=1:num_EEGchns;
            % Fast Fourier and Thompson's Multitaper Method Initializations
            nw=4; nfft_fft = 512;   nfft_pmtm = 2^nextpow2(size(oneTrial(:,chn),1));
            if nfft_pmtm < 100; % keep out extremely small trial sizes
                continue
            end
            
            % calculate the Fourier Transform (using Thompson's Multi-Taper Method)
                [BehaviorSegments(cl).Trials(p).pmtm_psd(:,chn),...
                    BehaviorSegments(cl).Trials(p).pmtm_confid(:,:,chn),...
                    BehaviorSegments(cl).Trials(p).pmtm_freq(:,chn)]...
                    = pmtm(oneTrial(:,chn), nw, nfft_pmtm, fs, 0.95);                
                
            % calculate the Fourier Transform (using Fast Fourier Algorithm)    
                BehaviorSegments(cl).Trials(p).fft(:,chn) = ...
                    fft(oneTrial(:,chn),nfft_fft)/size(oneTrial(:,chn),1);

            % calculate the short-time Fast Fourier transform (STFT)   
                wndw =20;     novrlp = wndw-1;      % initialize
                % Initialization Onset transitions    
                    Onset_trial = oneInitOnset; init_spect_flg=1;
                % Completion Onset transitions
%                     Onset_trial = oneCompOnset; comp_spect_flg=1;
                OnsetTriallength = length(Onset_trial);  % spectrogram number of frequency samples is the length of the trial period
                nfft_s = 2^nextpow2(OnsetTriallength);
                [ BehaviorSegments(cl).Trials(p).SPECTROGRAM(chn).stft,...
                    BehaviorSegments(cl).Trials(p).SPECTROGRAM(chn).freq, ...
                    BehaviorSegments(cl).Trials(p).SPECTROGRAM(chn).time,...
                    BehaviorSegments(cl).Trials(p).SPECTROGRAM(chn).stftpsd ] = ...
                    spectrogram(Onset_trial(:,chn), wndw, novrlp, nfft_s, fs); % generate STFT spectrogram
            end % repeat for all channels
        end     % repeat for all trials for one class
    end  % repeat for all six segmented classes of intended actions
    PREPROCESS.BehaviorSegments = BehaviorSegments;
    disp('____behavioral trials per class selected')

%% Channel Initializations
%{
    chnCPz = find(strcmp(PREPROCESS.EEG.channelOrder_rmvdCHNS,'CPz')); % to find the correct index for 'CPz'
    chnFC5 = find(strcmp(PREPROCESS.EEG.channelOrder_rmvdCHNS,'FC5')); % to find the correct index for 'C3'
    chnPOz = find(strcmp(PREPROCESS.EEG.channelOrder_rmvdCHNS,'POz')); % to find the correct index for 'POz'
    chnCP2 = find(strcmp(PREPROCESS.EEG.channelOrder_rmvdCHNS,'CP2')); % to find the correct index for 'CP2'
    selectedCHNS = [chnCPz, chnFC5, chnPOz];
%}
%% Set directory paths for 'topoplot_woDipole' function
%{
    addpath(genpath('\\172.27.216.40\Contreras-UH\Zach Hernandez\MATLAB\eeglab12_0_0_0b'))
    rmpath(genpath('\\172.27.216.40\Contreras-UH\Zach Hernandez\MATLAB\eeglab12_0_0_0b\functions\octavefunc'))
%}
%% Ranking Features by Channel and Frequency Range
%{
    freqbandList = [1 3; 3 6; 6 9; 9 12; 12 15; 15 18];
    num_freqbands = size(freqbandList,1);
    FeatureRanking = struct('GrandMean_AvgPwr',{},'GrandStdev_AvgPwr',{});
    for cl = 1:num_classes;
        for chn = [1:num_EEGchns];    
            for fb = 1:num_freqbands;
                num_trials = length(Behaviors(cl).Trials);
                
                % Initializing array of onset trial segments (per task-based class)
                avgPowerperTrial = zeros(1,num_trials);
                
                for p=1:num_trials;
                    oneTrial = BehaviorSegments(cl).Trials(p).EEGsignal;
                    nfft_pmtm = 2^nextpow2(size(oneTrial(:,chn),1));
                    if nfft_pmtm < 100; % keep out extremely small trial sizes
                        continue
                    end

                    % Initialize frequencies and PSDs for each trial   
                    oneTrial_freq = BehaviorSegments(cl).Trials(p).pmtm_freq(:,chn);
                    oneTrial_psd = BehaviorSegments(cl).Trials(p).pmtm_psd(:,chn);
                
                    % Initialize each frequency band to compute average power 
                    freqband_ind = find(freqbandList(fb,1) <= oneTrial_freq & ...
                        oneTrial_freq <= freqbandList(fb,2));
                    freqband_freq = oneTrial_freq(freqband_ind);
                
                    % Computing average power per trial per channel (using Thompson's Multi-Taper Method)
                    BehaviorSegments(cl).Trials(p).avgPower(chn) = trapz(freqband_freq,...
                        pow2db(oneTrial_psd(freqband_ind)));
                    avgPowerperTrial(p) = Behaviors(cl).Trials(p).avgPower;
                                   
                end     % repeat for all trials of one task-based class
                % Add the mean and standard deviation of calculated average
                % power for each channel and frequency group into array
                FeatureRanking(cl).GrandMean_AvgPwr(chn,fb) = mean(avgPowerperTrial);
                FeatureRanking(cl).GrandStdev_AvgPwr(chn,fb) = std(avgPowerperTrial);
            end     % repeat for all frequency bands  
        end         % repeat for all channels
    
        %Create a Topology Head plot for each frequency band
        keptchns = (1:64);          keptchns(chnsremoved) = [];  % generate indices of all channels kept for further analysis

        AvgPwr_64chns = zeros(64,1);
        for fb = (6); %1:num_freqbands;
            for chn = 1:num_EEGchns;
                AvgPwr_64chns(keptchns(chn)) = FeatureRanking(cl).GrandMean_AvgPwr(chn,fb);
            end
            %
            % Plot average power of each channel on a scalp topology map
            figure;    topoplot_woDipole((AvgPwr_64chns),chanLocs,'numcontour',0,...
            'plotchans',(keptchns),'electrodes','ptslabels','maplimits','maxmin',...
            'shading','interp','emarker',{'.','k',[],1},'gridscale',100,'whitebk','on');
            title(['Average Power Across ',tasklabels{cl},' Trials',10,'Frequency Group: ',...
                num2str(freqbandList(fb,1)),' - ',num2str(freqbandList(fb,2)),' Hz'],...
                'FontSize',16,'FontWeight','bold');
            xlim([-0.7 0.7]);   ylim([-0.7 0.7]);   set(gcf,'color','w');  
            caxis([-100 -15]);         CB_hndl = colorbar;      %For 6-9Hz caxis([-100 -15]);  
        	set(get(CB_hndl,'title'),'string',['Average',10,'Power (dB)'],'FontSize',16,'FontWeight','bold');
            set(CB_hndl,'FontSize',14);
            %set(CB_hndl,'Position',[0.948 0.0811 0.014 0.770]);
            %
        end         % repeat for all frequency bands
    end             % repeat for all task-based classes    
    % Take the ratio of both tasks
    FeatureRanking_Ratio = FeatureRanking(2).GrandMean_AvgPwr-FeatureRanking(1).GrandMean_AvgPwr;
    
   
    % Sort by Channels for each Frequency band
    [SortingbyChn,Index_chnsort] = sort(FeatureRanking_Ratio,'descend');
    % Plot values for each Frequency Band
    %
    for fb = (6);   %1:num_freqbands;
        
        AvgPwrRatio_64chns = zeros(64,1);
        for chn = 1:num_EEGchns;
                AvgPwrRatio_64chns(keptchns(chn)) = FeatureRanking_Ratio(chn,fb);
        end
        %
        % Plot average power of each channel on a scalp topology map
        figure;   topoplot_woDipole((AvgPwrRatio_64chns),chanLocs,'numcontour',0,...
            'plotchans',(keptchns),'electrodes','ptslabels','maplimits','maxmin',...
            'shading','interp','emarker',{'.','k',[],1},'gridscale',100,'whitebk','on');
        title(['Grand Mean Average Power Difference',10,'(Imitate-Observe)'],...
            'FontSize',16,'FontWeight','bold');
        ylabel(['Frequency Group: ',num2str(freqbandList(fb,1)),' - ',...
            num2str(freqbandList(fb,2)),' Hz'],'FontSize',14,'FontWeight','bold');   
        xlim([-0.7 0.7]);      ylim([-0.7 0.7]);    set(gcf,'color','w');
        CB_hndl = colorbar;
        set(get(CB_hndl,'title'),'string',['Average Power (dB)'],'FontSize',16,'FontWeight','bold');
        set(CB_hndl,'FontSize',14); set(CB_hndl,'Location','South');
        set(CB_hndl,'Position',[0.357 0.114 0.327 0.051]);
        %
        %Plot ratio values for each frequency group
        figure;  bar(SortingbyChn(:,fb)); 
        set(gca,'XTick',(1:num_EEGchns),'XTickLabel',PREPROCESS.EEG.channelOrder_rmvdCHNS(Index_chnsort(:,fb)),'FontSize',12,'FontWeight','bold');
        ylabel(['Average Power Difference',10,'(Imitate-Observe)'],'FontSize',14); xlim([1-0.5 num_EEGchns+0.5]);
        title(['Frequency Group: ',num2str(freqbandList(fb,1)),' - ',num2str(freqbandList(fb,2)),' Hz'],'FontSize',16);
    end
    %   
%}
%% Plotting for each Method of Analysis 
cd(InfantData_dir)
    %% Histograms of each class
        %   Only using CPz electrode as a sample signal for observation of
        %   histograms
        %{
        close all
    for chn = [chnCP2, chnFC5, chnPOz]; figure; % create a set of plots per channel
        set(gcf,'color','w');
        for cl = 1:num_classes;
            subplot(num_classes,1,cl);     num_trials = length(Behaviors(cl).Trials);

            Class_chnHist=zeros(length(histx),num_trials); 
            for p=1:num_trials; 
                Class_chnHist(:,p) = BehaviorSegments(cl).Trials(p).histogram(:,chn);
            end
            CLASShistogram = mean(Class_chnHist,2);
            bar(histx,CLASShistogram./max(CLASShistogram(:))); colormap(gray);
            
            text(2.95, 0.7, ['',num2str(sRound(kurt(CLASShistogram),3))],...
                'FontWeight','bold', 'FontSize',14, 'LineStyle','none','HorizontalAlignment','right');
        
            ylabel(tasklabels(cl),'FontSize',12,'FontWeight','bold');
            xlim([-3.2 3.2])
            set(specgraph.baseline,'LineWidth',3.0); 
            set(gca,'FontSize',16,'FontWeight','bold','LineWidth',2.0);
%             if cl==1; title(['Grand Mean Histograms across Session Trials: ',...
%                     PREPROCESS.EEG.channelOrder_rmvdCHNS{chn} ' electrode'],'FontSize',...
%                     20,'FontWeight','bold'); 
%             end
            if cl==num_classes; xlabel('Z-Score (std)','FontSize',14,'FontWeight','bold');
            else
                set(gca,'XTickLabel',[]);
            end
        end     %  repeat for all six task-based classes
    end   % repeat for all channels
    disp('____histograms constructed')
        %}
    %% Power Spectral Density of each class
        %   Only using CPz electrode as a sample signal for observation of
        %   power spectral density plots
        %{
        TaskcolorList = {'-r', '-b','-k','-g','-m','-c'};       mainLineHandleSet = zeros(1,num_classes);
        for chn = [chnCP2, chnFC5, chnPOz]; figure; % create a set of plots per channel     
            for cl = 1:num_classes;
%                 subplot(num_classes,1,cl);  
                
                num_trials = length(Behaviors(cl).Trials);

                % Initializing array of onset trial segments (per task-based class)
                f = fs/2*linspace(0,1,nfft_fft/2+1);
                ClassTrials_PSD = zeros(length(f),num_trials);
   
                for p=1:num_trials;
                    oneTrial = BehaviorSegments(cl).Trials(p).EEGsignal;
                    
                    %Concatenating the PSD of each trial into one array
                    ClassTrials_PSD(:,p) = 2*abs(BehaviorSegments(cl).Trials(p).fft(1:nfft_fft/2+1,chn));
                end     % repeat for all trials of one task-based class
                
                hndl_FFT(cl) = shadedErrorBar(f,ClassTrials_PSD',{@mean,@std},TaskcolorList{cl},1);  
                set(hndl_FFT(cl).mainLine,'LineWidth',2); xlim([0 14]); 
                ylim([min(mean(ClassTrials_PSD,2)-std(ClassTrials_PSD,0,2)) max(mean(ClassTrials_PSD,2)+std(ClassTrials_PSD,0,2))]);
                
                ylabel(['Spectral Power (\muV/Hz)'],'FontSize',14,'FontWeight','bold');
                
                if cl==1; title(['Grand Mean PSDs across Session Trials',10,...
                    PREPROCESS.EEG.channelOrder_rmvdCHNS{chn},' electrode'],...
                    'FontSize',16,'FontWeight','bold');
                end
                
                if cl==num_classes; xlabel('Frequency (Hz)','FontSize',14,'FontWeight','bold'); end
                mainLineHandleSet(cl) = hndl_FFT(cl).mainLine;                
                hold on
            end  %  repeat for all task-based classes
            hlgd = legend(mainLineHandleSet,tasklabels);
            set(hlgd,'FontSize',14,'FontWeight','bold','EdgeColor','w','Orientation','horizontal','Location','South')
            set(gca,'FontSize',12,'FontWeight','bold');
            set(gcf,'Color','w');
        end    % repeat for all channels
        
        %}
    %% Event-Related Potentials Along Initialization Onset of each class
        %   Plot the average of all EEG signal-based trials for each task
        %{
        for chn = [chnCP2, chnFC5, chnPOz]; figure; % create a set of plots per channel     
            for cl = 1:num_classes;
                subplot(num_classes,1,cl);  
                
                num_trials = length(Behaviors(cl).Trials);

                % Initializing time sample sizes
                trialTIME = linspace(-(OnsetTriallength/(2*fs)),(OnsetTriallength/(2*fs)),OnsetTriallength);

                % Initializing array of onset trial segments (per task-based class)
                Class_erp=zeros(length(Onset_trial),num_trials); 

                for p=1:num_trials; 
                    Class_erp(:,p) = BehaviorSegments(cl).Trials(p).initOnset(:,chn);
                end 
                hndl_OT = shadedErrorBar(trialTIME,Class_erp',{@mean,@std},'-r',1);  
                set(hndl_OT.mainLine,'LineWidth',2); xlim([-2 2]);  ylim([min(Class_erp(:)) max(Class_erp(:))]);
                hold on; line([0 0],[min(Class_erp(:)) max(Class_erp(:))], 'Color', 'k', 'LineWidth', 2,'LineStyle','--');
                
                ylabel([tasklabels(cl),num2str(num_trials),' Trials'],'FontSize',12,'FontWeight','bold');
                
                if cl==1; title(['Grand Mean ERPs across Session Trials',10,...
                    PREPROCESS.EEG.channelOrder_rmvdCHNS{chn},' electrode'],...
                    'FontSize',16,'FontWeight','bold');
                end
                
                if cl==num_classes; xlabel('Time (s)','FontSize',12,'FontWeight','bold'); end
            end  %  repeat for all six task-based classes
        end    % repeat for all channels
        %}
    %% Short-Time Fourier Transform Spectrogram of each class
        %   Only using CPz electrode as a sample signal for observation of
        %   power spectral density plots
        %{
        tic;
        for chn = [14 23]; figure; % create a set of plots per channel 
            for cl = 1:num_behaviors;
                subplot(num_behaviors,1,cl);            
                num_trials = length(Behaviors(cl).Trials);
%                 num_restTrials = size(PREPROCESS.CLASS.EEGclassTiming(6).EEGclassTrials,2);
            
                % Initialize time and frequency vectors
                spect_time = BehaviorSegments(cl).Trials(1).SPECTROGRAM(1).time-(OnsetTriallength/(2*fs)); % subtract out one second to center time at onset = zero seconds
                spect_freq = BehaviorSegments(cl).Trials(1).SPECTROGRAM(1).freq;

                %Initialize array of spectrogram trials to average
                Class_spect = zeros(length(spect_freq),length(spect_time), num_trials);
%                 Class_spectREST = zeros(length(spect_freq),length(spect_time), num_restTrials); 

                % concatenate all STFT-PSD trials for averaging
                for p=1:num_trials;
                    Class_spect(:,:,p) = db(BehaviorSegments(cl).Trials(p)...
                        .SPECTROGRAM(chn).stftpsd);
                end

                % concatenate all STFT trials of 'rest' for averaging
%                 for u=1:num_restTrials;
%                     Class_spectREST(:,:,u) = db(BehaviorSegments(cl).Trials(p)...
%                         .SPECTROGRAM(chn).stftpsd);
%                 end
                
                % Average all trials
                CLASS_spect_mean = mean(Class_spect,3);                
%                 CLASS_spectREST_mean = mean(Class_spectREST,3);
                [ALLTIME_PSD, ~] = pwelch(EEGsignalInput(:,chn),...
                    wndw*10, novrlp*10, nfft_s, fs);  
                
                % Subtract STFT-PSDs using 'rest' as a baseline
%                 CLASS_spect_mean_centered = CLASS_spect_mean-CLASS_spectREST_mean;
                CLASS_spect_mean_centered = CLASS_spect_mean - ...
                    db(repmat(ALLTIME_PSD,1,size(CLASS_spect_mean,2)));
                
                % Plot each spectrogram
                surf(spect_time, spect_freq, CLASS_spect_mean_centered,'LineStyle','none');
                view(0,90); axis tight; caxis([-20 20])
                Z_plotOverlay = max(CLASS_spect_mean_centered(:));
                hold on; line([median(spect_time) median(spect_time)],[0 20],[Z_plotOverlay Z_plotOverlay], 'Color', 'k', 'LineWidth', 2,'LineStyle','--');
                
                xlim([min(spect_time) max(spect_time)]);        ylim([0 20]);
                
                ylabel([tasklabels(cl),num2str(num_trials),' Trials'],'FontSize',12,'FontWeight','bold');
                if cl==1; title(['Grand Mean STFT Spectrograms across Session Trials, '...
                        ,num2str(wndw),' samples/window, ',num2str(round(100*(novrlp/wndw))),'% overlap',10,...
                        PREPROCESS.EEG.channelOrder_rmvdCHNS{chn},' electrode'],...
                        'FontSize',16,'FontWeight','bold'); 
                end
                
                if cl == num_behaviors-1 && exist('init_spect_flg','var'); 
                    xlabel('Initialization Onset Time (s)','FontSize',12,'FontWeight','bold');
                elseif cl == num_behaviors-1 && exist('comp_spect_flg','var'); 
                    xlabel('Completion Onset Time (s)','FontSize',12,'FontWeight','bold'); 
                end
            end  % repeat for all six segmented classes of intended actions
            CB_hndl = colorbar; set(get(CB_hndl,'title'),'string',['Spectral',10,' Power (dB)']);
            set(get(CB_hndl,'title'),'fontsize',14,'fontweight','bold');
            set(CB_hndl,'Position',[0.938 0.108 0.014 0.7614]);
        end % repeat for all three selected channels
        toc;
    %}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Generate FFT-PSD Spectrograms
        % Spectrogram Initializations
        %
            %             red       magenta      green       yellow         cyan       light purple
        taskIDColorList=[1 0 0;      1 0 1;      0 1 0;      1 1 0;         0 1 1;       1 0.8 1];
        tasklabels = {'Explore',   'Imitate',   'Point',  'Reach-Grasp', 'Reach-Offer',   'Rest'};
        wndw = 1024; novrlp = wndw-102; nfft=wndw; F = nfft*8; % (0.1:0.01:40);
        TimeFreqAnalysis=struct('FFTPSD',{},'COHERE',{},'COHERxy',{},'coherF',{},'CSPSD',{},'CPSDxy',{},'cspsdF',{});
        %}
    %{    
    % For the Three EEG Electrodes CPz, FC5, and POz %
        display('Generating PSD of FFT Spectrogram......');
        figure(20); SP_hndl = tight_subplot(4,1,[.01 .03],[.08 .01],[.05 .08]);
        for k=1:3;
            display(['     :: For ' num2str(PREPROCESS.EEG.channelOrder_rmvdCHNS{selectedCHNS(k)}) '......']);
            [~,FREQ,TIME,TimeFreqAnalysis(k).FFTPSD]...
                = spectrogram(transpose(PREPROCESS.EEG.resamp(:,selectedCHNS(k))),wndw,novrlp,F,fs); % for EEG
        % Plot the spectrogram       
            axes(SP_hndl(k)); surf(TIME,FREQ,db(TimeFreqAnalysis(k).FFTPSD),'LineStyle','none');   
            view(0,-90); axis tight; caxis([-80 80]); 
            text(-8,33,PREPROCESS.EEG.channelOrder_rmvdCHNS{selectedCHNS(k)},'FontWeight','bold');
        % Segment Trial into Classes using Vertical Lines
           for cl1 = 1:size(PREPROCESS.CLASS.TASKSEGMENTS.data,1);
                X_start = min(TIME(abs(TIME - PREPROCESS.CLASS.TASKSEGMENTS.data(cl1,1)/100) <= diff(TIME(1:2))/2));
                X_end = min(TIME(abs(TIME - PREPROCESS.CLASS.TASKSEGMENTS.data(cl1,2)/100) <= diff(TIME(1:2))/2));
                Z_plotOverlay = min(db(TimeFreqAnalysis(k).FFTPSD(:)));
                colorSelect = taskIDColorList(PREPROCESS.CLASS.TASKSEGMENTS.data(cl1,3),:);

                line([X_start X_start], [0.1 40], [Z_plotOverlay Z_plotOverlay], 'Color', 'k', 'LineWidth', 1.5);
                line([X_start-0.125 X_start+0.125], [20 20], [Z_plotOverlay Z_plotOverlay],...
                    'Color', colorSelect, 'LineWidth', 3);
                line([X_end X_end], [0.1 40], [Z_plotOverlay Z_plotOverlay], 'Color', 'k', 'LineWidth', 1.5);
                line([X_end-0.25 X_end+0.25], [20 20], [Z_plotOverlay Z_plotOverlay],...
                    'Color', colorSelect, 'LineWidth', 3);
            end
            %            
        end
        ylabel('Frequency (Hz)','FontSize',14);  
        CB_hndl = colorbar; set(get(CB_hndl,'title'),'string',['Spectral',10,' Power (dB)']);
        set(CB_hndl,'Position',[0.948 0.350 0.014 0.557]);
%--------------------------------------------------------------------------------------------------------------------    
    % For the Head Acceleration Sensor %
        k=4;
        display('     :: For Head Accel Sensor......');
        [~,~,~,TimeFreqAnalysis(k).FFTPSD]...
            = spectrogram(transpose(PREPROCESS.KINE.InfantForehead.resampGCMA),wndw,novrlp,F,fs); % for Head Accel Sensor            
    % Plot the spectrogram
        figure(20);  axes(SP_hndl(k));
        surf(TIME,FREQ,db(TimeFreqAnalysis(k).FFTPSD),'LineStyle','none'); view(0,-90); axis tight;
        caxis([-50 50]); 
        text(-10,33,['HEAD' 10 'ACCEL'],'FontWeight','bold');           
    % Segment Trial into Classes using Vertical Lines
        line_hdl1=zeros(size(PREPROCESS.CLASS.TASKSEGMENTS.data,1),1);
        task_hdl1=cell(size(PREPROCESS.CLASS.TASKSEGMENTS.data,1),1);
        for cl1 = 1:size(PREPROCESS.CLASS.TASKSEGMENTS.data,1);
            X_start = min(TIME(abs(TIME - PREPROCESS.CLASS.TASKSEGMENTS.data(cl1,1)/100) <= diff(TIME(1:2))/2));
            X_end = min(TIME(abs(TIME - PREPROCESS.CLASS.TASKSEGMENTS.data(cl1,2)/100) <= diff(TIME(1:2))/2));
            Z_plotOverlay = min(db(TimeFreqAnalysis(k).FFTPSD(:)));
            task_hdl1{cl1} = tasklabels{PREPROCESS.CLASS.TASKSEGMENTS.data(cl1,3)};
            colorSelect = taskIDColorList(PREPROCESS.CLASS.TASKSEGMENTS.data(cl1,3),:);
            
            line([X_start X_start], [0.1 40], [Z_plotOverlay Z_plotOverlay], 'Color', 'k', 'LineWidth', 1.5);
            line_hdl1(cl1) = line([X_start-0.125 X_start+0.125], [20 20], [Z_plotOverlay Z_plotOverlay],...
                'Color', colorSelect, 'LineWidth', 3);
            line([X_end X_end], [0.1 40], [Z_plotOverlay Z_plotOverlay], 'Color', 'k', 'LineWidth', 1.5);
            line([X_end-0.25 X_end+0.25], [20 20], [Z_plotOverlay Z_plotOverlay],...
                'Color', colorSelect, 'LineWidth', 3);
        end
        MarkerColorSet = get(line_hdl1, 'color');
        MarkerColorSet=cell2mat(MarkerColorSet);
        [uniqueMarkerColorSet,ia,~] = unique(MarkerColorSet,'rows');
        legend(line_hdl1(ia),task_hdl1(ia),'Orientation', 'Horizontal','Location',...
            'SouthOutside','Color',[0.5 0.5 0.5],'EdgeColor',[0.8 0.8 0.8],'Position',[0.0156 0.0145 0.363 0.0293]);

        CB_hndl = colorbar; set(get(CB_hndl,'title'),'string',['Spectral',10,' Power (dB)']); 
        set(CB_hndl,'Position',[0.948 0.0811 0.014 0.18]);
        xlabel('Time (s)','FontSize',14);
        set(SP_hndl,'XLim',[675 775],'YLim',[0.1 40],'YScale','log','XTick',(690:20:770),'YTick',[0.2,4,40],...
            'YTickLabel',[0.2,4,40],'YDir','normal','FontWeight','bold','FontSize',12);
        set(SP_hndl(1:3),'XTickLabel',{});
    %}
    %% Generate Wavelet Transform Coherence Time-Frequency Plots (by Grinsted)    
        
    %% Generate Squared-Magnitude Coherence Spectrogram (Zach's Code)
    %{
    display('Generating Magnitude Coherence Spectrogram......');
    addpath(strcat(pwd,'\ShortTimeCoherenceSpectrogram'))
    y = PREPROCESS.KINE.InfantH.resampGCMA; 
        figure(2); SP_hndl = tight_subplot(3,1,[.01 .03],[.08 .01],[.05 .08]);
    for k=1:3;
        display(['     :: For Head Accel Sensor to ' num2str(PREPROCESS.EEG.channelOrder_rmvdCHNS{selectedCHNS(k)}) '......']);
        x = PREPROCESS.EEG.resamp(:,selectedCHNS(k));
        [TimeFreqAnalysis(k).COHERE,FREQ1,TIME1] = spectrogram_mscohere(x,y,wndw,novrlp,nfft,fs); 
    % Plot the spectrogram
    %
        axes(SP_hndl(k)); surf(TIME1,FREQ1,TimeFreqAnalysis(k).COHERE,'LineStyle','none'); view(0,-90); axis tight; 
        text(max(TIME1)+5,max(FREQ1)-5,['H.A. + ' PREPROCESS.EEG.channelOrder_rmvdCHNS{selectedCHNS(k)}],'FontWeight','bold');
        caxis([0 1]); CB_hndl = colorbar;
        if k==1; set(get(CB_hndl,'title'),'string',['Squared Magnitude',10,'Coherence (A.U.)']);
            set(CB_hndl,'Position',[0.948 0.0811 0.014 0.770]);   end
        if k==2; ylabel('Frequency (Hz)','FontSize',14);end        
        if k==3; xlabel('Time (s)','FontSize',14); end
        
    % Segment Trial into Classes using Vertical Lines
        for cl1 = 1:size(PREPROCESS.CLASS.TASKSEGMENTS.data,1);
            X_start = min(TIME1(abs(TIME1 - PREPROCESS.CLASS.TASKSEGMENTS.data(cl1,1)/100) <= diff(TIME1(1:2))/2));
            X_end = min(TIME1(abs(TIME1 - PREPROCESS.CLASS.TASKSEGMENTS.data(cl1,2)/100) <= diff(TIME1(1:2))/2));
            Z_plotOverlay = max(db(TimeFreqAnalysis(k).COHERE(:)));
            colorSelect = taskIDColorList(PREPROCESS.CLASS.TASKSEGMENTS.data(cl1,3),:);
            
            line([X_start X_start], [0.1 40], [Z_plotOverlay Z_plotOverlay], 'Color', 'k', 'LineWidth', 1.5);
            line([X_start-0.125 X_start+0.125], [20 20], [Z_plotOverlay Z_plotOverlay],...
                'Color', colorSelect, 'LineWidth', 3);
            line([X_end X_end], [0.1 40], [Z_plotOverlay Z_plotOverlay], 'Color', 'k', 'LineWidth', 1.5);
            line([X_end-0.25 X_end+0.25], [20 20], [Z_plotOverlay Z_plotOverlay],...
                'Color', colorSelect, 'LineWidth', 3);
        end
         line([675 775],[1 1],[Z_plotOverlay Z_plotOverlay], 'Color', 'w', ...
             'LineStyle','--', 'LineWidth',3);
         line([675 775],[4 4],[Z_plotOverlay Z_plotOverlay], 'Color', 'w', ...
             'LineStyle','--', 'LineWidth',3);
    end
    set(SP_hndl,'XLim',[675 775],'YLim',[0.1 40],'YScale','log','XTick',(690:20:770),'YTick',[0.2,4,40],...
        'YTickLabel',[0.2,4,40],'YDir','normal','FontWeight','bold','FontSize',12);
    set(SP_hndl(1:2),'XTickLabel',{}); 
    %}

    %% Generate Squared-Magnitude Coherence Spectrogram (Atilla's Code) 
    %{
    w_sizeSpectgram=1024;                        w_sizeCohere=w_sizeSpectgram/8;    shift=1;    
    noverlapSpectgram=w_sizeSpectgram-shift;    noverlapCohere=w_sizeCohere-1;      nfft=2^nextpow2(w_sizeCohere);
    sig1=OPAL(trial).accel_resamp(1:end-abs(length(EEG(trial).signal_resamp)-length(OPAL(trial).accel_resamp)),:); 
    ntotaltime=size(sig1,1);
    total_nwindows = fix((ntotaltime-noverlapSpectgram)/(w_sizeSpectgram-noverlapSpectgram));
    tic; figure(2); SP_hndl = tight_subplot(3,1,[.01 .03],[.08 .01],[.05 .08]);
    for k=1:3;
        sig2=EEG(trial).signal_filtered_zscore(:,selectedCHNS(k));     
%         if length(F)<1; nfft = w_sizeSpectgram; F=nfft/2+1; end
        cohxy=zeros(length(F),total_nwindows);      FREQ1=zeros(length(F),total_nwindows);

        id1=1;
        temp1=zeros(w_sizeSpectgram,total_nwindows);
        temp2=zeros(w_sizeSpectgram,total_nwindows);
        for i=1:total_nwindows;
            temp1(:,i)=sig1(id1:id1+w_sizeSpectgram-1,1); 
            temp2(:,i)=sig2(id1:id1+w_sizeSpectgram-1,1);     
            id1=id1+shift;
        end
%         wb = waitbar(0,'Atilla`s Coherence spectrogram being constructed. Please wait...');
        progressStepSize = 100; ppm = ParforProgMon('Calculating coherence spectrogram array: ', total_nwindows, progressStepSize, 300, 80);
        parfor i=1:total_nwindows;
            display(['Iteration ' num2str(i)])
            [cohxy(:,i),FREQ1(:,i)] = mscohere_noPERMISSIONS(temp1(:,i),temp2(:,i),w_sizeCohere,noverlapCohere,F,fs);
            if mod(i,progressStepSize)==0; ppm.increment(); end
%             waitbar(i/100)
        end
%         close(wb)
        ppm.delete()
        TimeFreqAnalysis(k).COHERE = cohxy;
        TIME1 = (((0:(total_nwindows-1))*(w_sizeSpectgram-noverlapSpectgram))+((w_sizeSpectgram)/2)')/fs; 
    % Plot the spectrogram      
        axes(SP_hndl(k)); surf(TIME1,FREQ1,TimeFreqAnalysis(k).COHERE,'LineStyle','none'); view(0,-90); axis tight; 
        text(max(TIME1)+5,max(FREQ1)-5,['H.A. + ' channels{selectedCHNS(k)}],'FontWeight','bold');
        caxis([0 1]); CB_hndl = colorbar;
        if k==1; set(get(CB_hndl,'title'),'string','Coherence (A.U.)'); end
        if k==2; ylabel('Frequency (Hz)','FontSize',14);end        
        if k==3; xlabel('Time (s)','FontSize',14); end
    % Segment Trial into Classes using Vertical Lines
        for cl1 = 2:num_classes;
            line([TIME1(abs(TIME1 - EEG(trial).signal_classes(cl1).time(1)) <= diff(TIME1(1:2))/2)...
                TIME1(abs(TIME1 - EEG(trial).signal_classes(cl1).time(1)) <= diff(TIME1(1:2))/2)],...
                [min(FREQ1) max(FREQ1)],...
                [min(TimeFreqAnalysis(k).COHERE(:)) min(TimeFreqAnalysis(k).COHERE(:))],'Color','k','LineWidth',3);
            line([TIME1(abs(TIME1 - EEG(trial).signal_classes(cl1).time(end)) <= diff(TIME1(1:2))/2)...
                TIME1(abs(TIME1 - EEG(trial).signal_classes(cl1).time(end)) <= diff(TIME1(1:2))/2)],...
                [min(FREQ1) max(FREQ1)],...    
                [min(TimeFreqAnalysis(k).COHERE(:)) min(TimeFreqAnalysis(k).COHERE(:))],'Color','k','LineWidth',3);           
        end
    end
    %}
%     display(['Spectrogram computation time:: ' num2str(toc/60) ' minutes or ' num2str(toc) ' seconds'])
    %% Generate Cross Spectrum PSD Spectrogram
    %{
    display('Generating Cross Spectrum PSD Spectrogram......');
    y = OPAL(trial).accel_resamp(1:end-1,:);
    figure(3); SP_hndl = tight_subplot(3,1,[.01 .03],[.08 .01],[.05 .08]);
    for k=1:3;
        x = EEG(trial).signal_resamp(:,selectedCHNS(k));
        [TimeFreqAnalysis(k).CSPSD,FREQ1,TIME1] = spectrogram_cpsd(x,y,wndw,novrlp,F,fs); %F only used for spectrogram windows, nfft=wndw used for cross PSD calculations
    % Plot the spectrogram
        magnCSPSD = db(abs(TimeFreqAnalysis(k).CSPSD));     anglCSPSD = -angle(TimeFreqAnalysis(k).CSPSD); 
        axes(SP_hndl(k)); surf(TIME1,FREQ1,magnCSPSD,'LineStyle','none'); view(0,-90); axis tight; 
        text(max(TIME1)+5,max(FREQ1)-5,['H.A. + ' channels{selectedCHNS(k)}],'FontWeight','bold');
        caxis([80 -80]); CB_hndl = colorbar;
        if k==1; set(get(CB_hndl,'title'),'string',['Cross Spectrum',10,'Phase (rad)']); end
        if k==2; ylabel('Frequency (Hz)','FontSize',14);end        
        if k==3; xlabel('Time (s)','FontSize',14); end
    % Segment Trial into Classes using Vertical Lines
        for cl1 = 2:num_classes;
            line([TIME1(abs(TIME1 - EEG(trial).signal_classes(cl1).time(1)) <= diff(TIME1(1:2))/2)...
                TIME1(abs(TIME1 - EEG(trial).signal_classes(cl1).time(1)) <= diff(TIME1(1:2))/2)],...
                [min(FREQ1) max(FREQ1)],...
                [min(db(abs(TimeFreqAnalysis(k).CSPSD))) min(db(abs(TimeFreqAnalysis(k).CSPSD)))],'Color','k','LineWidth',3);
            line([TIME1(abs(TIME1 - EEG(trial).signal_classes(cl1).time(end)) <= diff(TIME1(1:2))/2)...
                TIME1(abs(TIME1 - EEG(trial).signal_classes(cl1).time(end)) <= diff(TIME1(1:2))/2)],...
                [min(FREQ1) max(FREQ1)],...    
                [min(db(abs(TimeFreqAnalysis(k).CSPSD))) min(db(abs(TimeFreqAnalysis(k).CSPSD)))],'Color','k','LineWidth',3);           
        end
    end
    set(SP_hndl,'YLim',[min(FREQ1) max(FREQ1)],'YScale','log','YTick',[0.2,4,40],...
        'YTickLabel',[0.2,4,40],'YDir','normal');
    set(SP_hndl(1:2),'XTickLabel',{}); 
    %}  
        %display('All spectrograms for all channels completed');        
    %% Save the Spectrogram
    %{
        display('Saving Spectrogram');
        dataStorage_dir=['\\bmi-nas-01\Contreras-UH\Zach Hernandez\MATLAB\MAIN_LFDA_GMM\EEGcoherence_Results\S' num2str(subject)];
        mkdir(dataStorage_dir)
        cd(dataStorage_dir)
%         Figure_Name_tif=strcat('S',num2str(subject),'T',num2str(trial),'_raw3CHNandHEADSENSOR_CoherenceSpectrogram.tif');
%         saveas(gca,Figure_Name_tif);
        Figure_Name_fig=strcat('S',num2str(subject),'T',num2str(trial),'_raw3CHNandHEADSENSOR_CoherenceSpectrogram.fig');
        saveas(gca,Figure_Name_fig);
        display('Spectrogram saved');
    %}
%-------------------------------------------------------------------------%
%% Display Computation Time
display(['TOTAL TIME:: ' num2str(toc/60) ' minutes or ' num2str(toc) ' seconds'])