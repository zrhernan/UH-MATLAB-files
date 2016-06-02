%% for calling helper functions
addpath('\\172.27.216.40\Contreras-UH\Infantdata\Infantdata\Source-Estimation\uh_fieldtrip')
addpath('\\172.27.216.40\Contreras-UH\Infantdata\Infantdata\Source-Estimation\uh_fieldtrip\demos')
addpath('\\172.27.216.40\Contreras-UH\Infantdata\Infantdata\Source-Estimation\uh_fieldtrip\fileio')
addpath('\\172.27.216.40\Contreras-UH\Infantdata\Infantdata\Source-Estimation\uh_fieldtrip\preprocess')
addpath('\\172.27.216.40\Contreras-UH\Infantdata\Infantdata\Source-Estimation\uh_fieldtrip\time_freq')
addpath('\\172.27.216.40\Contreras-UH\Infantdata\Infantdata\Source-Estimation\uh_fieldtrip\utils')
addpath('\\172.27.216.40\Contreras-UH\Infantdata\Infantdata\Source-Estimation\uh_fieldtrip\visualization')
addpath('\\bmi-nas-01\Contreras-UH\Infantdata\Infantdata\code\Zachs_Infant_decoding_files\WaveletCoherenceToolbox-grinsted')
addpath('\\bmi-nas-01\Contreras-UH\Infantdata\Infantdata\code\Zachs_Infant_decoding_files')
addpath(genpath('C:\Users\zrhernan\Documents\MATLAB\fieldtrip'))

%% Generate List of Infant Data Folders
InfantDir = '\\172.27.216.40\Contreras-UH\Infantdata\Infantdata\Data\';
disp_flag = 0;
[ InfantDataAnnotList, InfantID ] = defineInfantFolders( InfantDir, disp_flag );

%%
for ii=18%:length(InfantID) %starting from CM17
    disp(['Opening EEG recording of subject ',InfantID{ii}])
    infantfolder = InfantDataAnnotList{ii};
    serverPath1 = [InfantDir,infantfolder];
    cd(serverPath1)
    
    %% Load EEG FT dataset (w/trials)
    % Check if the EEG file exists
    fullFileName1 = 'EEG\FORCeCleanedEEG_noWindowing_removedhighZchns.mat';
    if ~exist(fullFileName1, 'file')
      % File does not exist.
      warningMessage = sprintf('Warning: file does not exist:\n%s', fullFileName1);
      disp(warningMessage)
      disp('Skipping to next infant data set')
      continue
    else
%         load('EEG\importedRawEEG.mat','data')
%         load('EEG\rejectedchans_resampled_preprocEEG.mat','dataRejectHighZchn');
%         load('EEG\EEGTriggerTimes.mat');
%         load('EEG\ExperStartStopIndices.mat');
%         shift1 = eventTriggers(ExperStartInd);
%         data.time{1} = data.time{1}(1:end-shift1+1);
%         data.trial{1} = data.trial{1}(:,shift1:end);
%         data.sampleinfo(2) = length(data.time{1});
        load(fullFileName1) % for FORCe method
        eegset = dataClean; clear dataClean;
%         ASRdata=load('EEG\rejectedchans_resampled_preprocEEG_ASR.mat'); % for ASR method
%         eegset = dataRejectHighZchn; clear dataRejectHighZchn;
        eegset.sampleinfo = eegset.cfg.previous.trl(:,1:2);
        eegset.cfg.conditionlabel = eegset.cfg.previous.conditionlabel;
    end
    
    %% Load Kinematics FT dataset (w/trials)   
    % Check if the kinematics file exists
    fullFileName2 = 'Kinematics\resampled500HzKINEtrl.mat';
    if ~exist(fullFileName2, 'file')
      % File does not exist.
      warningMessage = sprintf('Warning: file does not exist:\n%s', fullFileName2);
      disp(warningMessage)
      disp('Skipping to next infant data set')
      continue
    else
%         load('Kinematics\importedRawOPALdata_FTstructure.mat','datakine_filt')
%         load('Kinematics\ExperStartStopIndices.mat');
%         load('Kinematics\KINETriggerTimes.mat');
%         shift2 = TrigEdgeOnset_1000Hz(ExperStartInd);
%         datakine_filt.time{1} = datakine_filt.time{1}(1:end-shift2+1);
%         datakine_filt.trial{1} = datakine_filt.trial{1}(:,shift2:end);
%         datakine_filt.sampleinfo(2) = length(datakine_filt.time{1});
        load(fullFileName2,'datakineTrlresamp')
        accelset = datakineTrlresamp; clear datakineTrlresamp;
    end
    
    %% Divide up trials by class
    eeg_dClass = uh_classextract(eegset);
    kine_dClass = uh_classextract(accelset);

    %% Wavelet coherence (WTC)
    % The WTC finds regions in time frequency space where the two
    % time series co-vary (but does not necessarily have high power).
    
    % for whole session coherence (warning: takes too long)
    %{
    fs = 1000;       % both are sampled at 1000 Hz
    ts1 = data.time{1}; % extract time samples that correspond to the data
    ts2 = datakine_filt.time{1}; % extract time samples that correspond to the data
    tsmin = min(length(ts1),length(ts2)); % minimum number of time sample length between each time series
    ts = data.time{1}(1:tsmin)'; %common time sampling    
    chncntr = 0;
    for chn = 1:length(data.label)
        channame = data.label{chn};
        if any(strcmp(channame,{'Fz'}))==1,%'CPz','Oz'
            chncntr = chncntr + 1;
            acccntr = 0;
%             for acc = 7%:9%:length(datakine_filt.label)
%                 accname = datakine_filt.label{acc};
%                 acccntr = acccntr + 1;
                
                ds1 = transpose(zscore(data.trial{1}(chn,1:tsmin)));
                ds2 = zscore(sqrt(sum(datakine_filt.trial{1}(7:9,1:tsmin)'.^2,2)));
                options = {'mcc',10};
                tic;
                [Rsq,period,~,coi,wtcsig,Wxy] = wtc([ts ds1],[ts ds2],options{:}); % to speed up Monte Carlo test
                timecompute(toc);
%             end
            figure(chncntr), set(gcf,'color',[1 1 1])
            subplot(3,1,acccntr)
            makeWTCfigure(Rsq,period,coi,wtcsig,Wxy,ts)
%             title(['WTC: ' channame ' & ' accname] )
        else
            continue
        end
    end
    %}
    
    % Using trials divided by class
    %{
    fs = 500;       % both are sampled at 500 Hz
    ts = eeg_dClass{1}.time{1}; % extract time samples that correspond to the data
    clscntr = 0;
    for cls = 3%:length(eeg_dClass)
        clscntr = clscntr + 1;
        classname = eeg_dClass{cls}.conditionlabel;
        chncntr = 0;
        for chn = 1:length(eeg_dClass{cls}.label)
            channame = eeg_dClass{cls}.label{chn};
            acccntr = 0;
            if any(strcmp(channame,{'Fz','CPz','Oz'}))==1
                chncntr = chncntr + 1;
                for acc = 7:9%:size(kine_dClass{cls}.trial{trl},1)
                    accname = kine_dClass{cls}.label{acc};
                    acccntr = acccntr + 1;
                    Rsq = cell(1,length(eeg_dClass{cls}.trial));
                    coi = cell(1,length(Rsq));
                    wtcsig = cell(1,length(Rsq));
                    Wxy = cell(1,length(Rsq));
                    for trl = 120%:length(eeg_dClass{cls}.trial)
                        trialname = ['Trial ',num2str(trl)];
                        ds1 = zscore(eeg_dClass{cls}.trial{trl}(chn,:));
                        ds2 = zscore(kine_dClass{cls}.trial{trl}(acc,:));
                        options = {'mcc',100};
                        [Rsq{trl},period,~,coi{trl},wtcsig{trl},Wxy{trl}] = wtcwPARFOR([ts;ds1]',[ts;ds2]',options{:}); % to speed up Monte Carlo test
                    end
                    figure(clscntr), set(gcf,'color',[1 1 1])
                    subplot(3,3,sub2ind([3,3],acccntr,chncntr))
                    MeanRsq = mean(cat(3,Rsq{:}),3);
                    MeanWxy = mean(cat(3,Wxy{:}),3);
                    Meanwtcsig = mean(cat(3,wtcsig{1:20}),3);
                    makeWTCfigure(MeanRsq,period,coi{1},Meanwtcsig,MeanWxy,ts)
                    makeWTCfigure(Rsq{trl},period,coi{trl},wtcsig{trl},Wxy{trl},ts)
                    title(['WTC: ' classname ',' trialname ': ' channame ' & ' accname] )
                end
            else
                continue
            end
        end
    end
%}
    
    % Using trials not divided by class
%{
    fs = 500;       % both are sampled at 500 Hz
    ts = eegset.time{1}; % extract time samples that correspond to the data
    chncntr = 0;
    for chn = 1:length(eegset.label)
        channame = eegset.label{chn};
        acccntr = 0;
        if any(strcmp(channame,{'Fz','CPz','Oz'}))==1
            chncntr = chncntr + 1;           
            Rsq = cell(1,length(eegset.trial));
            coi = cell(1,length(Rsq));
            wtcsig = cell(1,length(Rsq));
            Wxy = cell(1,length(Rsq));
            for trl = 1:20%length(eegset.trial) %120
                trialname = ['Trial ',num2str(trl)];
                ds1 = zscore(eegset.trial{trl}(chn,:));
                ds2 = zscore(sqrt(sum(accelset.trial{trl}(7:9,:)'.^2, 2)))';
                options = {'mcc',100};
                [Rsq{trl},period,~,coi,wtcsig{trl},Wxy{trl}] = wtcwPARFOR([ts;ds1]',[ts;ds2]',options{:}); % to speed up Monte Carlo test
            end
            
            MeanRsq = mean(cat(3,Rsq{:}),3);
            MeanWxy = mean(cat(3,Wxy{:}),3);
            Meanwtcsig = mean(cat(3,wtcsig{1:20}),3);
            
            figure(1), set(gcf,'color',[1 1 1])
            subplot(3,1,chncntr)
            makeWTCfigure(MeanRsq,period,coi{1},Meanwtcsig,MeanWxy,ts)
            makeWTCfigure(Rsq{trl},period,coi{trl},wtcsig{trl},Wxy{trl},ts)
            title(['WTC of ' trialname ': ' channame ' & Head Magnitude Acceleration'] )
        else
            continue
        end
    end
%}

    % Using signals averaged across all trials not divided by class
%{
    fs = 500;       % both are sampled at 500 Hz
    ts = eegset.time{1}; % extract time samples that correspond to the data
    MagnAccelset = zeros(length(accelset.trial),length(ts));
    headsensidcs = find(strncmp(accelset.label,'InfantsForehead',15)==1);
    for trl = 1:length(accelset.trial),
        MagnAccelset(trl,:) = sqrt(sum(accelset.trial{trl}(headsensidcs,:)'.^2, 2));
    end
    MuMagnAccelset = mean(MagnAccelset);
    MuEEGset = mean(cat(3,eegset.trial{:}),3);
    SelectedChannelList = {'Fz','CPz','Oz'};
    chncntr = 0;
    for chn = 1:length(eegset.label)
        channame = eegset.label{chn};
        if any(strcmp(channame,SelectedChannelList))==1
            chncntr = chncntr + 1;
            ds1 = zscore(MuEEGset(chn,:));
            ds2 = zscore(MuMagnAccelset);
            options = {'mcc',100};
%             if strcmp(InfantID{ii},'CM17'), options = {options,'AR1','burg' };
            [Rsq,period,~,coi,wtcsig,Wxy] = wtcwPARFOR([ts;ds1]',[ts;ds2]',options{:}); % to speed up Monte Carlo test
            
            figure(1), set(gcf,'color',[1 1 1])
            subplot(3,1,find(ismember(SelectedChannelList,channame)==1))
            makeWTCfigure(Rsq,period,coi,wtcsig,Wxy,ts)
            title(['WTC of ' channame ' & Head Magnitude Acceleration'] )
        else
            continue
        end
    end
%}
    % Using signals averaged across all trials and divided by class
%
    fs = 500;       % both are sampled at 500 Hz
    ts = eegset.time{1}; % extract time samples that correspond to the data
    headsensidcs = find(strncmp(accelset.label,'InfantsForehead',15)==1);
    SelectedChannelList = {'Fz','CPz','Oz'};

    clscntr = 0;
    for cls = 1:length(eeg_dClass)
        clscntr = clscntr + 1;
        classname = eeg_dClass{cls}.conditionlabel;
        chncntr = 0;
        MagnAccelset = zeros(length(kine_dClass{cls}.trial),length(ts));
        for trl = 1:length(kine_dClass{cls}.trial),
            MagnAccelset(trl,:) = sqrt(sum(kine_dClass{cls}.trial{trl}(headsensidcs,:)'.^2, 2));
        end
        MuMagnAccelset = mean(MagnAccelset);
        MuEEGset = mean(cat(3,eeg_dClass{cls}.trial{:}),3);
        if all(isempty(MuEEGset)) || all(isnan(MuMagnAccelset)), continue, end
        ts = eeg_dClass{cls}.time{1}; % extract time samples that correspond to the data
        for chn = 1:length(eegset.label)
            channame = eegset.label{chn};
            if any(strcmp(channame,SelectedChannelList))==1
                chncntr = chncntr + 1;
                ds1 = zscore(MuEEGset(chn,:));
                ds2 = zscore(MuMagnAccelset);
                options = {'mcc',100,'S0',1/40,'ms',1/.001};
    %             if strcmp(InfantID{ii},'CM17'), options = {options,'AR1','burg' };
                [Rsq,period,~,coi,wtcsig,Wxy] = wtcwPARFOR([ts;ds1]',[ts;ds2]',options{:}); % to speed up Monte Carlo test

                figure(1), set(gcf,'color',[1 1 1])
                subplot(6,3,sub2ind([3,6],find(ismember(SelectedChannelList,channame)==1),clscntr))
                makeWTCfigure(Rsq,period,coi,wtcsig,Wxy,ts,options)
                if clscntr == 1, title([channame ' & |ACCEL_H_E_A_D|']), end
%                 if strcmp('Oz',channame), colorbar, end
                if strcmp('Fz',channame), ylabel([upper(classname),10,10,'Frequency (Hz)']), end
                if clscntr == 6, xlabel('Time (sec)'), end
                set(gca,'FontSize',12,'FontWeight','bold')
            else
                continue
            end
        end
    end
    cbhndl = colorbar;
    title(cbhndl,['Magnitude', 10,'Coherence'],'fontsize',14)
    set(cbhndl,'position',[0.935,0.111,0.02,0.818])
%}
%% Save it!
myfig=figure(1);
enhanceFigure(myfig)
filename = 'wtcohere_meanEEGchans_headaccel.tif';
cd([serverPath1,'\Kinematics'])
print(myfig,'-dpng', filename, '-r600')
close all
end