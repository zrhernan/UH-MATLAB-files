  % 
  % matlab file InfantDataForwardSelectionMethod
  % 
  % Author1: Jesus G Cruz-Garza (jesusc90@gmail.com) 
  % Author2: Zachery R Hernandez
  % Date:    April 2014
  % 
  % Function   : Forward Selection of EEG Channels that contribute the 
  % most to classification. 
  % 
  % Description: Identify and sort the EEG channels that contributed the
  % most to the classification of the Infant Data using the LFDA-GMM
  % algorithm.
  % 
  % 
  % 
  % 
  % 






clc; close all; clear all; format compact
tic  % start computation time

%% Annotation for Infant Data
InfantDataAnnotList = {'N09-12-04-2013','J20-08-19-2013','B06-10-30-2013',...
    'GR09-07-12-2014','A06-09-28-2013','LW10-06-19-2014'};
infant = menu(['Which infant would you',10,' like to analyze?'],...
    'N09','J20','B06','GR09','A06','LW10');
close

%% Initialization of Data Structures
PROCESS = struct('EEG',{},'KINE',{},'CLASS',{},'CLASSIFIER',{});
disp('Initialization of Data Structures')

%% Initializing and Assigning directory paths....
InfantData_dir = 'C:\Users\zrhernan\Infant_decoding_files\';
cd(InfantData_dir)
disp('Initializing and Assigning directory paths')
toc

%% Importation of Class Labeling Vector
disp('Importing Class Labeling Vector')
load(['Data\',InfantDataAnnotList{infant},'\classlabel_7-classes.mat']);     % load class vector
% load(['Data\',InfantDataAnnotList{infant},'\classlabel_Observe&Imitate.mat']);     % load class vector
PROCESS(1).CLASS.classlabel=classlabel; clear classlabel % and save into data structure
disp('____Target vector imported')
toc

%% Importation of List of Time-Segmented Tasks
LastColumnIndex = [104,147,68,114,49];              LastColumnIndex_2class = [39,46,19]; %for binary class dataset
LastColumnIndex_7C_J20 = 162; 
TASKSEGMENTS = importTaskTrialInfo(['Data\',InfantDataAnnotList{infant},'\7-class start-stop times.txt'],...
    2, LastColumnIndex_7C_J20); %LastColumnIndex(infant));  %    % load list of start and stop times for each class
% TASKSEGMENTS = importTaskTrialInfo(['Data\',InfantDataAnnotList{infant},'\Observe&Imitate_start-stop times.txt'],...
%     2, LastColumnIndex_2class);  % load list of start and stop times for each 'observe' and 'imitate' class
PROCESS.CLASS.TASKSEGMENTS = TASKSEGMENTS;   % and save into PROCESS structure
disp('____List of tasks imported')
toc

%% Importation of EEG Data
PROCESS.EEG = load(['Data\',InfantDataAnnotList{infant},'\EEGfiles.mat']);  % load structure of EEG attributes

disp('____EEG data imported')

toc
%% Resampling Data to 100 Hz
% resampling for EEG data {originally sampled at 1000 Hz}      
num_EEGchns = size(PROCESS.EEG.raw_vidsync,2);
for k=1:num_EEGchns;
    PROCESS.EEG.resamp(:,k) = decimate(PROCESS.EEG.raw_vidsync(:,k),10); 
end

% Initializing time sample sizes
fs=100; %[Hz] sampling rate of the decimated signal
tmax=size(PROCESS.EEG.resamp,1)/fs; % number of time samples
t=0:1/fs:tmax; t1=transpose(t(1:end-1)); % EEG time vector
disp(['____ALL data resampled to ' num2str(fs) ' Hz'])  
toc

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
    delta_freq=[1 4]; % band passing into the delta frequency
    PROCESS.EEG.filtered_cell = filter_data_bpass(PROCESS.EEG.resamp_cell, fs, n_f, delta_freq);
    PROCESS.EEG.filtered = cell2mat(PROCESS.EEG.filtered_cell);
    disp('____EEG data filtered')
    toc
    
%%  Lag-Based Feature Extraction
    % for a sampling frequency of  100 Hz, the lags will vary from 0 to -90 ms in steps of 10 ms
    lag = (0:9); max_lag=10; %maximum lag of -90 ms 
    PROCESS.EEG.featmat = timelag(PROCESS.EEG.filtered,lag);     % generate feature matrix
    PROCESS.CLASS.classlabel=PROCESS.CLASS.classlabel(max_lag:end);     % truncate target vector to maintain same number of rows
    disp('____Feature matrix produced')
    toc
    
%%  Standardization across channels (Z-scores: subtracting mean from each data point and deviding by standard deviation )
    %tip: mean 'should be already' zero after filtering, check this if your filter is working correctly
    %PROCESS.EEG.zscores_featmat = zscore(PROCESS.EEG.featmat);
    for i = 1:size(PROCESS.EEG.featmat,2)
        mn = mean(PROCESS.EEG.featmat(:,i)); ss = std(PROCESS.EEG.featmat(:,i));
        PROCESS.EEG.zscores_featmat(:,i) =  (PROCESS.EEG.featmat(:,i) - mn)/ss;
    end
    disp('____Feature matrix standardized')
    toc
    
%% LFDA-GMM Clasifier Algorithm

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

%========================================================================            
            

disp('initializing FORWARD SELECTION ALGORITHM')
%%
%=======================================================================
%======================FORWARD SELECTION ALGORITHM======================
%=======================================================================
% ======================================================================
% Description
% ======================================================================
% A forward selection algorithm. It runs the classifier (using a fixed set
% of 10 dimensions) for each EEG channel, and selects the channel that
% yields the highest classification accuracy. The selected channel is
% paired with every other channel at a time. The best paired channel is
% selected. The algorithm continues until all channels are selected. 
%
% ======================================================================
% Input / Output If we look at this section as a function
% ======================================================================
% Input: 
%       tr_perc         % array of training percentage values to use
%       dimSet          % set of values for optimization of LDFA  
%                       parameter 'dim'
%       knnSet          % set of values for optimization of LDFA 
%                       parameter 'knn'
%       maxiterations   % maximum number of iterations to run for the 
%                       LFDA GMM classifier due to random sampling. 
%                       More than 10 iterations preferred
%       num_EEGchns     % number of EEG channels considered
%
% Output: 
%       SelectedList_MeanSTD    % mean (for 'maxiterations' number of 
%                               iterations) standard deviation of the
%                               output of the classifier with the channels
%                               selected.
%       SelectedList_MeanACC    % mean classification accuracy of the
%                               output of the classifier with the channels
%                               selected
%       SelectedChansList       % list, in order, of channels selected
%                               that yield the highest classification
%                               accuracy with the input values specified
%
% ======================================================================
% Other Variables
% ======================================================================
% SelectedList_MeanSTD: This list will be updated after every iteration
% with the standard deviation of the channel that produced the highest 
% accuracy. 
%
% SelectedList_MeanACC=[]: This list will be updated after every iteration
% with the mean accuracy of the channel that produced the highest accuracy.
%
% SelectedChnsList=[]: This list will be updated after every iteration
% with the channel number that produced the highest accuracy.
%
% SelectedChnFeatures = []: This list is used to select the columns of the 
% feature matrix that correspond to all the channels selected.
%
% dimSet: Fixed number of dimensions to have as output for the 
% dimensionality reduction (LFDA)technique. 
%
% knnSet: Fixed number of k-nearest neighbors.
%
% maxiterations: number of interations to run the LFDA-GMM algorithm for
% each EEG channel. The mean and standard deviation of the accuracy is then
% calculated.
%
% EEGchns_leftover_count = 1; 
%
%
% num_EEGchns: maximum number of EEG channels.
%
% fwsind: Forward selection indice. It goes from 1 to the number of
% channels being analyzed. 
%
% EEGchns_leftover: Channels that have not yet been selected. 
%
% mean_accuracy: vector of mean classification accuracy yielded when each
% individual channel is tested paired with SelectedChansList. 
%
% std_accuracy: vector of standard deviation of the classification 
% accuracies yielded when each individual channel is tested paired with
% SelectedChansList.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Jesus G. Cruz-Garza, University of Houston, September 2014.
FWD_SELECTION = struct('SelectedChnFeatures',{},'SelectedChnsList',{},...
    'SelectedList_MeanACC',{},'SelectedList_StdACC',{});

FWD_SELECTION(1).SelectedChnFeatures = {};
FWD_SELECTION.SelectedChnsList = [];
dim = dimSet;
knn = knnSet;

for FWDSelection_step = 1:num_EEGchns
    EEGleftoverchns_indices = 1:num_EEGchns;
% Remove channels already selected in previous foward selection iterations 
    chnsList = [];
    for chns = 1:FWDSelection_step - 1
        if FWDSelection_step == 1
            chnsList = [];
        else
            chnsList(chns) = FWD_SELECTION(chns).SelectedChnsList;
        end
    end
    EEGleftoverchns_indices(chnsList) = [];
    
    disp('=======================');
    disp(' --------------------');
    disp('Selected List of Channels')
    disp(chnsList')
    disp('=======================');
    EEGchns_leftover_count = 1;
    
    for selectedchn_index = EEGleftoverchns_indices
        FWD_SELECTION(FWDSelection_step).SelectedChnFeatures{selectedchn_index} = [];
        
        disp(['+++Progress+++: testing channel ', num2str(selectedchn_index),...
            ' with Selected Channels List: ',...
            num2str(chnsList')])

                    %%====================
                    %%=======LFDA GMM=====
                    %%====================

        % Selecting indices of the full feature matrix for classification
        length_SCL = length(FWD_SELECTION); % length of forward selection steps taken thus far
        for i = 1:length_SCL
            firstFeature = FWD_SELECTION(i).SelectedChnsList*10-9;
            lastFeature = FWD_SELECTION(i).SelectedChnsList*10;
            FWD_SELECTION(FWDSelection_step).SelectedChnFeatures{selectedchn_index} = ...
                [FWD_SELECTION(FWDSelection_step).SelectedChnFeatures{selectedchn_index},...
                (firstFeature:lastFeature)];
        end
        FWD_SELECTION(FWDSelection_step).SelectedChnFeatures{selectedchn_index} = ...
            [FWD_SELECTION(FWDSelection_step).SelectedChnFeatures{selectedchn_index},...
            (selectedchn_index*10-9:selectedchn_index*10)];

        currentFeatIndices = FWD_SELECTION(FWDSelection_step).SelectedChnFeatures{selectedchn_index};

        disp('feature indices used ')
        disp(currentFeatIndices);
        feat_size = size(currentFeatIndices,2);
        disp([num2str(feat_size),' features selected']);

        % Redefine new feature matrix and target class vector using
        % selected indices
        ForwardSelectionFeatMat = FeatMatXsubj(:,currentFeatIndices);
        ForwardSelectionClassLbl = ClassLblXsubj;
    
        % Initializations
        classmatrixlist = cell(1,maxiterations);
        confusionmatrixlist = cell(1,maxiterations);
        accuracylist = zeros(maxiterations,1);
        precisionlist = cell(1,maxiterations);
        sensitivitylist = cell(1,maxiterations);
        
        % open Parallel Computing Toolbox
        if matlabpool('SIZE') == 0;   matlabpool('open');   end
        
        parfor iteration = 1:maxiterations
            disp(['...Progress: Foward Selection Step = ',...
                num2str(FWDSelection_step), ', chns = ',...
                num2str(EEGchns_leftover_count), '/', ...
                num2str((num_EEGchns+1)-FWDSelection_step),...
                '  ',num2str(iteration), ' iterations'])

            disp(['...tested channel ', num2str(selectedchn_index),...
                ' w/ Selected Channels List = ',...
                num2str(chnsList')])

        %% Random Sub-Sampling Cross Validation
            trainprct=tr_perc;    
        % Providing a percentage of training and testing samples to use for the classifier
            [data1,C1] = ytoc(ForwardSelectionFeatMat',ClassLblXsubj);
        % percentage of training data samples to use (using percentage of least populated class as reference)    
            num_train=round(min(C1)*trainprct/100); 
        % number of testing data samples to use (smallest class size - number of training samples)
            num_test=min(C1)-num_train;
        % randomly split dataset into samples set for training and testing 
            [data_train,data_test,c_train,c_test] = datasplitfcn(data1((C1(1):end),:),C1(2:end),'2',num_train,num_test);
        % assign training data and corresponding class labels per time sample
            X = data_train';
            Y = ctoy(c_train);

        %% Local Fisher's Discriminant Analysis (LFDA)
            % dim=12; knn=1 ; %knn=odd number     
            display('LFDA')
                                                         %    dim:        dimensionality of reduced space (default: d)
            [T,Z] = lfda(X,Y,dim,'plain',knn);           %    T  : d x r  transformation matrix (Z=T'*X)
                                                         %    Z  : r x n  matrix of dimensionality reduced samples 
                                                         %              metric: type of metric in the embedding space (default: 'weighted')
                                                         %              'weighted'        --- weighted eigenvectors 
                                                         %              'orthonormalized' --- orthonormalized
                                                         %              'plain'           --- raw eigenvectors

            if ~isreal(T);
                disp('Error: LDFA resulted in complex eigenvalues')
                continue
            end
        %% Gaussian Mixture Modeling (GMM) 
            ext_data_train = Z;     % training data for GMM
            ext_data_test  = data_test * T;     % testing data for GMM
            max_cluster = 10;  % Maximum number of clusters considered by GMM
            reg = 1e-5;        % Regularization parameter used in GMM
            display('GMM')
            try
                [class_matrix, density, posterior, comp_weight, obj, value] = gmmclassifier(ext_data_train,  ext_data_test,  c_train,  c_test,  max_cluster,  reg);
            catch err
                disp('Error within the GMM classifier');
                disp(err.getReport)
                continue 
            end
            classmatrixlist{iteration}=class_matrix;
            confusion_matrix = confusionmatrix_withSpecificities(class_matrix);      % compute confusion matrix with specificity per class included
            confusionmatrixlist{iteration} = confusion_matrix;
            accuracylist(iteration) = confusion_matrix(end,end-1);        % vector used to save accuracies
            precisionlist{iteration} = confusion_matrix(end,1:end-2);
            sensitivitylist{iteration} = confusion_matrix(1:end-1,end-1);

            %update2=['Iteration number ',num2str(iteration),',and accuracy: ',num2str(accuracylist(iteration))];

        end % repeat for each iteration of the LFDA-GMM classifier (for random sub-sampling X-validation)
        % matlabpool('close') % close Parallel Computing Toolbox
        % add LFDA-GMM outputs to data structure 'PROCESS'
        PROCESS.CLASSIFIER(dimIDX,knnIDX).classmatrixlist=classmatrixlist;
        PROCESS.CLASSIFIER(dimIDX,knnIDX).confusionmatrixlist = confusionmatrixlist;
        PROCESS.CLASSIFIER(dimIDX,knnIDX).accuracylist = accuracylist;       
        PROCESS.CLASSIFIER(dimIDX,knnIDX).precisionlist=precisionlist;
        PROCESS.CLASSIFIER(dimIDX,knnIDX).sensitivitylist=sensitivitylist;

        FWD_SELECTION(FWDSelection_step).mean_accuracy(selectedchn_index) = mean(accuracylist);
        FWD_SELECTION(FWDSelection_step).std_accuracy(selectedchn_index) = std(accuracylist);

        disp('accuracy update')
        disp(FWD_SELECTION(FWDSelection_step).mean_accuracy(selectedchn_index))
        EEGchns_leftover_count = EEGchns_leftover_count + 1;
        %} 
    end % repeat for all channels leftover (not in the classification accuracy-based selected list)
    %% End of LFDA GMM
    
    %% Select  best performing channel and add it to the 'selected channels list'
%     all_meanAccuracies = diag(FWD_SELECTION(FWDSelection_step).mean_accuracy);
%     all_meanStandardDeviations = diag(FWD_SELECTION(FWDSelection_step).std_accuracy);
    
    [accsort, indice] = sort(FWD_SELECTION(FWDSelection_step).mean_accuracy,'descend');
    
    FWD_SELECTION(FWDSelection_step).SelectedList_MeanACC = accsort(1);
    
    FWD_SELECTION(FWDSelection_step).SelectedList_StdACC = FWD_SELECTION(FWDSelection_step).std_accuracy(indice(1));
    
    FWD_SELECTION(FWDSelection_step).SelectedChnsList = indice(1);
    
    
    disp('====================');
    disp('++++++++++++++++++++');
    for steps = 1: length(FWD_SELECTION)
        disp(['Selected Channels: ',num2str(FWD_SELECTION(steps).SelectedChnsList)])
        disp(['Selected Accuracies: ',num2str(FWD_SELECTION(steps).SelectedList_MeanACC)])
    end
    %disp(['Selected Standard deviations:',num2str(FWD_SELECTION(steps).SelectedChnsList]_StdACC)])
    disp('++++++++++++++++++++');    
    disp('====================');
    disp('++++++++++++++++++++');    
    disp('====================');
end
matlabpool('close')
PROCESS.FWD_SELECTION = FWD_SELECTION; % add to processed data structure
%% Displaying
disp('Forward Selection Algorithm results:: ');
disp('Selected:    Channel:   Acc:');
for steps = 1: length(FWD_SELECTION)
disp(horzcat(steps,...
    FWD_SELECTION(steps).SelectedChnsList,...
    FWD_SELECTION(steps).SelectedList_MeanACC));
end
    

