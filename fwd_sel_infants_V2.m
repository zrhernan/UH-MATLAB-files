clc; close all; clear all; format compact

%% PRE-PROCESSING
%% Initialization of Data Structures
PROCESS = struct('EEG',{},'KINE',{},'CLASS',{},'CLASSIFIER',{});
disp('Initialization of Data Structures')

%% Initializing and Assigning directory paths....
InfantDataJohn_dir = '\\bmi-nas-01\Contreras-UH\Infantdata\Data\John-08-19-2013\';
LFDAGMMData_dir = '\\bmi-nas-01\Contreras-UH\Zach Hernandez\MATLAB\MAIN_LFDA_GMM\';
cd(InfantDataJohn_dir)
disp('Initializing and Assigning directory paths')

%% Importation of Class Labeling Vector
load('dataSegment\classlabel.mat');     % load class vector
PROCESS(1).CLASS.classlabel=classlabel; clear classlabel % and save into data structure
TASKSEGMENTS = importdata('video\class start-stop times.txt');  % load list of start and stop times for each class
PROCESS.CLASS.TASKSEGMENTS = TASKSEGMENTS;   % and save into PROCESS structure
disp('Importing Class Labeling Vector')
disp('____Target vector imported')

%% Importation of EEG Data
load ('EEG\john-08-19-2013.mat');  % load structure of EEG attributes
PROCESS(1).EEG = EEG; clear EEG % and save into data structure
PROCESS(1).EEG.raw_uV = double((PROCESS.EEG.data').*(PROCESS.EEG.gain));  % save EEG data into field of 'EEG' structure
for chn=1:length(PROCESS.EEG.chanlocs); PROCESS.EEG.channel_order(chn) = cellstr(PROCESS.EEG.chanlocs(chn).labels); end % save a cell list of channel names
% Missing from EEG set: 'FT9','FT10','PO9','PO10'  
% create list of electrodes with very high impedances  
highZchnlist = {'C3','Cz','P8','AF7','C5','C2','CP4','P2','P6','PO7','PO3','PO4','PO8'}; % Not in channel set: 'PO10'
% create list of electrodes along the periphery of the head
peripheralchnlist = {'Fp1','Fp2','AF8','F7','F8','FT7','FT8','T7','T8','TP9',...
    'TP7','TP8','TP10','P7','O1','Oz','O2'}; % Already in 'highZchnlist': 'P8','AF7','PO7','PO8'
removechnlist = {highZchnlist{:}, peripheralchnlist{:}};  % combine these two lists together
chnsremoved=zeros(length(removechnlist),1);
for rmv = 1:length(removechnlist); chnsremoved(rmv)=find(strcmp(PROCESS.EEG.channel_order,removechnlist{rmv})==1); end
PROCESS.EEG.raw_uV(:,chnsremoved)=[];   % and remove them from the EEG dataset
PROCESS.EEG.channel_order(:,chnsremoved)=[];   % and remove them from the EEG channels list
% Shifting Data to match timing of video recording
EEGStart=PROCESS.EEG.event(5).latency; % sychronized start of EEG data to video recording
EEGEnd=PROCESS.EEG.event(43).latency; % sychronized end of EEG data to video recording
PROCESS.EEG.raw_vidsync = PROCESS.EEG.raw_uV(EEGStart:EEGEnd,:); % new EEG set to reflect video recording times
% Acquiring Non-EEG data
PROCESS.EOG.raw = PROCESS.EEG.raw_uV(:,strcmp(PROCESS.EEG.channel_order,'EOG'));
PROCESS.ECG.raw = PROCESS.EEG.raw_uV(:,strcmp(PROCESS.EEG.channel_order,'ECG1')|strcmp(PROCESS.EEG.channel_order,'ECG2'));
% ... and remove from EEG array
nonEEGchannels = [find(strcmp(PROCESS.EEG.channel_order,'EOG')),find(strcmp(PROCESS.EEG.channel_order,'ECG1')),find(strcmp(PROCESS.EEG.channel_order,'ECG2'))];
PROCESS.EEG.raw_vidsync(:,nonEEGchannels)=[];   % and remove them from the EEG dataset
PROCESS.EEG.channel_order(:,nonEEGchannels)=[];   % and remove them from the EEG channels list
disp('____EEG data imported') 


    %% Initialize New Class of 'Don't Cares' between classes
    %{
        transition_width = 1000;
        classlabel_original=PROCESS.CLASS.classlabel;
        [PROCESS.CLASS.classlabel,classlabel_original] = IDCClass(PROCESS.CLASS.classlabel,...
            classlabel_original,transition_width);
        disp('____class transition samples removed')
    %}

%% Resampling Data to 100 Hz
    % resampling for EEG data {originally sampled at 1000 Hz}      
        num_EEGchns = size(PROCESS.EEG.raw_vidsync,2);
        for k=1:num_EEGchns;
            PROCESS.EEG.resamp(:,k) = decimate(PROCESS.EEG.raw_vidsync(:,k),10); 
        end
    % Precise alignment of both data streams
    % Take out extra samples in the beginning. Should only be 1-5 samples to remove
    samplesrmvd = abs(length(PROCESS.EEG.resamp) - length(PROCESS.CLASS.classlabel));
    if samplesrmvd > 5;
        disp('Resampled biomedical data streams larger than five samples, re-check syncing alignment.');
    else
        if length(PROCESS.EEG.resamp) > length(PROCESS.CLASS.classlabel);
            PROCESS.EEG.resamp = PROCESS.EEG.resamp(samplesrmvd+1:end,:);
        elseif length(PROCESS.CLASS.classlabel) > length(PROCESS.EEG.resamp);
            PROCESS.CLASS.classlabel = PROCESS.CLASS.classlabel(samplesrmvd+1:end,:);
        end
    end

    % Initializing time sample sizes
        fs=100; %[Hz] sampling rate of the decimated signal
        tmax=size(PROCESS.EEG.resamp,1)/fs; % number of time samples
        t=0:1/fs:tmax; t1=transpose(t(1:end-1)); % EEG time vector
        disp(['____ALL data resampled to ' num2str(fs) ' Hz'])   
%% Transform signal from MAT to CELL
%
% Initialization of cell array
    PROCESS.EEG.resamp_cell=cell(1,num_EEGchns);
    % For creating the cell array of EEG data
    for k=1:num_EEGchns;
        PROCESS.EEG.resamp_cell{1,k}=PROCESS.EEG.resamp(:,k);
    end
%%  Extract Delta Band EEG using a band-pass filter
    n_f=3;  % filter order
    delta_freq=[0.2 4]; % band passing into the delta frequency
    PROCESS.EEG.filtered_cell = filter_data_bpass(PROCESS.EEG.resamp_cell, fs, n_f, delta_freq);
    PROCESS.EEG.filtered = cell2mat(PROCESS.EEG.filtered_cell);
    disp('____EEG data filtered')
%%  Lag-Based Feature Extraction
    % for a sampling frequency of  100 Hz, the lags will vary from 0 to -90 ms in steps of 10 ms
    lag = (0:9); max_lag=10; %maximum lag of -90 ms 
    PROCESS.EEG.featmat = timelag(PROCESS.EEG.filtered,lag);     % generate feature matrix
    PROCESS.CLASS.classlabel=PROCESS.CLASS.classlabel(max_lag:end);     % truncate target vector to maintain same number of rows
    disp(['____Feature matrix produced']);
%%  Standardization across channels (Z-scores: subtracting mean from each data point and deviding by standard deviation )
    %tip: mean 'should be already' zero after filtering, check this if your filter is working correctly
    PROCESS.EEG.zscores_featmat = zscore(PROCESS.EEG.featmat);
    disp('____Feature matrix standardized')
%% LFDA-GMM Clasifier Algorithm
%     PROCESS(1).CLASSIFY = struct('mean_accuracy',{},'std_accuracy',{},'mean_precision',{},...
%         'std_precision',{},'mean_sensitivity',{},'std_sensitivity',{},'classmatrixlist',{});
%---Classifier Input Definitions-----------------------------------------------%
FeatMatXsubj=PROCESS.EEG.zscores_featmat;
ClassLblXsubj=PROCESS.CLASS.classlabel;
%------------------------------------------------------------------------------%
%---Optimization Initializations-----------------------------------------------%
    dimIDX=1;  % for indexing the dimension parameter
    knnIDX=1;  % for indexing the kNN parameter
    tpIDX=1;   % for indexing the percentage of training samples
    tr_perc = (50); % array of training percentage values to use
    dimSet = 10;  % set of values for optimization of LDFA parameter 'dim'
    knnSet = 7;      % set of values for optimization of LDFA parameter 'knn'
    maxiterations = 30; % maximum number of iterations to run for the LFDA GMM classifier due to random sampling. More than 10 iterations preferred 
%------------------------------------------------------------------------------%
% clear dat1 EEGdata dat chf
%========================================================================            
            



%%
%=======================================================================
%======================FORWARD SELECTION ALGORITHM======================
%=======================================================================
dim=dimSet; knn=knnSet;
sss=[];
A=[];
IX=[];
list = [];

for jjjo = 1:size(PROCESS.EEG.filtered,2)
    kept = 1:size(PROCESS.EEG.filtered,2);
    kept(IX) = [];
    mean_accuracy = zeros(1,length(kept));
    std_accuracy = zeros(1,length(kept));
    disp('=======================');
    disp(num2str(jjjo));
    disp('--------------------');
    disp(IX')
    disp('=======================');
    for onlythis = kept

        

            %% ===================
            %%=======LFDA GMM=====
            %%====================
 
           
            
for i = 1:length(IX)
a = IX(i)*10-9;
b= IX(i)*10;
list = [list a:b];
end
list = [list onlythis*10-9:onlythis];
size(IX)
size(list)

fmrx=FeatMatXsubj(:,list);
ic=ClassLblXsubj;


    trainprct=tr_perc;
        accuracylist=[];   %(length(XvalidList(:,1)),1);
        precisionlist=zeros(maxiterations,max(ic)-1); sensitivitylist=zeros(max(ic)-1,maxiterations);
        classmatrixlist=[];
        
            parfor (iteration = 1:maxiterations,30)   %length(XvalidList(:,1));
                
                  disp(iteration)

                % Providing a percentage of training and testing samples to
                  %     use for the classifier
                  [data1,C1] = ytoc(fmrx',ic);
                  td1=round(min(C1)*trainprct/100); %td2=round(min(C2)*trainprct/100); ed2=round(min(C2)*(100-trainprct)/100);
                  ed1=min(C1)-td1;
 
                  [data_train,data_test,c_train,c_test] = datasplitfcn(data1(1:end-C1(length(C1)),:),C1(1:end-1),'2',td1,ed1);

                    X = data_train';
                    Y = ctoy(c_train);

                %% LFDA-GMM        
                display('LFDA')
                                                             %    dim:        dimensionality of reduced space (default: d)
                [T,Z] = lfda(X,Y,dim,'plain',knn);           %    T  : d x r  transformation matrix (Z=T'*X)
                                                             %    Z  : r x n  matrix of dimensionality reduced samples 
                                                             %              metric: type of metric in the embedding space (default: 'weighted')
                                                             %              'weighted'        --- weighted eigenvectors 
                                                             %              'orthonormalized' --- orthonormalized
                                                             %              'plain'           --- raw eigenvectors

                ext_data_train = Z;
                ext_data_test  = data_test * T;
    %             ext_data_test = data2 * T;
    %             c_test=C2;
                max_cluster = 10;  % Maximum number of clusters considered by GMM
                reg = 1e-5;        % Regularization parameter used in GMM
                display('GMM')
                try
                [class_matrix, density, posterior, comp_weight, obj, value] = gmmclassifier(ext_data_train,  ext_data_test,  c_train,  c_test,  max_cluster,  reg);
                catch err
                    continue 
                end
                classmatrixlist=horzcat(classmatrixlist,class_matrix);
                confusion_matrix = confusionmatrixeditedZH(class_matrix);

                accuracylist(iteration) = confusion_matrix(end,end-1);        %vector used to save accuracies

                update2=['Iteration number ',num2str(iteration),',and accuracy: ',num2str(accuracylist(iteration))];

            end


            mean_accuracy(onlythis) = mean(accuracylist);
            std_accuracy(onlythis) = std(accuracylist);

            SummedConfusMat=confusionmatrixeditedZH(classmatrixlist);


disp('accuracy update')
disp(mean_accuracy)
disp('standard deviation accuracy update')
disp(std_accuracy)
disp('total iterations')
disp(maxiterations)




    
    disp('Forward Selection % Completion for current channel')
    disp(onlythis/37 *100)

        
    end
    %% End of LFDA GMM
    
    %% Select  best performing channel and add it to the 'selected channels list'
    [accsort, indice] = sort(mean_accuracy,'descend');
    sss= [sss;std_accuracy(indice(1))];
    A = [A;accsort(1)];
    IX = [IX;indice(1)];
    disp('====================');
    disp('++++++++++++++++++++');
    disp(A)
    disp(IX)
    disp(sss)
    disp('++++++++++++++++++++');    
    disp('====================');
    disp('FORWARD SELECTION % COMPLETION')
    disp(onlythis/37 *100)
    disp('++++++++++++++++++++');    
    disp('====================');
end


%% Displaying

disp('Forward Selection Algorithm results:: ')
disp(IX)

    
    

