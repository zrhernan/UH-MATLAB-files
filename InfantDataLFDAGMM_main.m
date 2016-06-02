clc, clear all, close all
tic;  % start computation time

% for calling helper functions
addpath(genpath('\\bmi-nas-01\Contreras-UH\Infantdata\Infantdata\code\Zachs_Infant_decoding_files'))

%% Annotation for Infant Data
InfantDataAnnotList = {'N09-12-04-2013','J20-08-19-2013','B06-10-30-2013',...
    'GR09-07-12-2014','A06-09-28-2013','LW10-06-19-2014','AR16-07-14-2014',...
    'RB23-12-04-2014','A18-08-15-2013'};

%% Selection of infants
% infant = menu(['Which infant would you',10,' like to analyze?'],...
%     'N09','J20','B06','GR09','A06','LW10','AR16','RB23');  %,'A18'
% close
InfantID = {'N09','J20','B06','GR09','A06','LW10','AR16','RB23','A18'};

for infant = 2% 1:length(InfantID);
%% Initialization of Data Structures
PROCESS = struct('EEG',{},'KINE',{},'CLASS',{},'FEATSELECT',{},...
    'CLASSIFIER',{},'LFDA_OPTIMIZATION',{},'VALIDATION',{});
% PROCESS(1).KINE = struct('InfantH',{},'InfantArmL',{},'InfantArmR',{},'ActorArmL',{},'ActorArmR',{});

%% Initializing and Assigning directory paths....
serverPath1 = ['\\bmi-nas-01\Contreras-UH\Infantdata\Infantdata\Data\',...
    InfantDataAnnotList{infant}];
serverPath2 = ['\\bmi-nas-01\Contreras-UH\Infantdata\Infantdata\code\',...
    'Zachs_Infant_decoding_files\Data\',InfantDataAnnotList{infant}];
desktopPath = ['C:\Users\zrhernan\Infant_decoding_files\Data\',...
    InfantDataAnnotList{infant}];
mkdir(serverPath1,'Classification Results')
cd([serverPath1 '/Classification Results'])
diary('ON');  % open diary file used for saving classification info
disp(['Performing Classification using EEG data from ', InfantID{infant}])
diary('OFF');  % close diary file used for saving classification info
cd(serverPath1)

%% Importation of Class Labeling Vector
load([serverPath1,'\Behavioral Segmentation\classlabel.mat']);     % load class vector
% load(['Data\',InfantDataAnnotList{infant},'\classlabel_7-classes.mat']);     % load class vector
% load(['Data\',InfantDataAnnotList{infant},'\classlabel_Observe&Imitate.mat']);     % load class vector
num_classes = max(classlabel); % define the numbers of behaviors used for classification
PROCESS(1).CLASS.classlabel=classlabel; clear classlabel % and save into data structure
disp('____Target vector imported')

%% Importation of List of Time-Segmented Task
%{
LastColumnIndex = [104,144,68,114,51,49,153,238,83];      LastColumnIndex_2class = [39,46,19]; %for binary class dataset
LastColumnIndex_7C_J20 = 162; 
TASKSEGMENTS = importTaskTrialInfo(['Data\',InfantDataAnnotList{infant},'\class start-stop times.txt'],...
    2, LastColumnIndex(infant)); %LastColumnIndex_7C_J20); %     % load list of start and stop times for each class
% TASKSEGMENTS = importTaskTrialInfo(['Data\',InfantDataAnnotList{infant},'\Observe&Imitate_start-stop times.txt'],...
%     2, LastColumnIndex_2class);  % load list of start and stop times for each 'observe' and 'imitate' class
PROCESS.CLASS.TASKSEGMENTS = TASKSEGMENTS;   % and save into PROCESS structure
disp('____List of tasks imported')
%}
%% Importation of EEG Data
PROCESS.EEG = load([serverPath1,'\EEG\EEGfiles.mat']);  % load structure of EEG attributes
disp('____EEG data imported')
disp('--------------------------------------------------')
timecompute(toc)
disp('--------------------------------------------------')

%% Importation of Kinematics Data
%{
PROCESS.KINE.InfantH = load (['Data\',InfantDataAnnotList{infant},'\InfantsForehead_SI-000773.mat']); % save kinematics Data into each body sensor
PROCESS.KINE.InfantT = load (['Data\',InfantDataAnnotList{infant},'\InfantsTrunk_SI-000742.mat']); % save kinematics Data into each body sensor
PROCESS.KINE.InfantLA = load (['Data\',InfantDataAnnotList{infant},'\InfantsLeftArm_SI-000722.mat']); % save kinematics Data into each body sensor
PROCESS.KINE.InfantRA = load (['Data\',InfantDataAnnotList{infant},'\InfantsRightArm_SI-000738.mat']); % save kinematics Data into each body sensor
PROCESS.KINE.FullSet = horzcat(PROCESS.KINE.InfantH.Acc_n_syncd,...
    PROCESS.KINE.InfantT.Acc_n_syncd,PROCESS.KINE.InfantLA.Acc_n_syncd,PROCESS.KINE.InfantH.Acc_n_syncd);
disp('____Kinematics data imported')
%}

%% Resampling Data to 100 Hz
    % resampling for EEG data {originally sampled at 1000 Hz}      
        num_EEGchns = size(PROCESS.EEG.uVdata_rmvdCHNS_syncd,2);
        for k = 1:num_EEGchns;
            PROCESS.EEG.resamp(:,k) = decimate(PROCESS.EEG.uVdata_rmvdCHNS_syncd(:,k),10); 
        end
%{
    % resampling for gravity-compensated magnitude acceleration (GCMA) data {originally sampled at 128 Hz}
        num_KINEsnsrs = size(PROCESS.KINE.FullSet_syncd,2);
        for p = 1:num_KINEsnsrs;     
            PROCESS.KINE.resamp = resample(PROCESS.KINE.FullSet_syncd,100,128); 
        end
    
    % Precise alignment of both data streams
    % Take out extra samples in the beginning. Should only be 1-5 samples to remove
%     samplesrmvd = abs(length(PROCESS.EEG.resamp) - length(PROCESS.KINE.resamp));
    if samplesrmvd > 5;
        disp('Resampled data streams larger than five samples, re-check syncing alignment.');
    else
        if length(PROCESS.EEG.resamp) > length(PROCESS.KINE.resamp);
            PROCESS.EEG.resamp = PROCESS.EEG.resamp(samplesrmvd+1:end,:);
        elseif length(PROCESS.KINE.resamp) > length(PROCESS.EEG.resamp);
            PROCESS.KINE.resamp = PROCESS.KINE.resamp(samplesrmvd+1:end,:);
        end
    end
%}
    % Initializing time sample sizes
        fs=PROCESS.EEG.srate/10; %[Hz] sampling rate of the decimated signal
        tmax=size(PROCESS.EEG.resamp,1)/fs; % number of time samples
        t=0:1/fs:tmax; t1=transpose(t(1:end-1)); % EEG time vector
%         t2 = linspace(0,(length(PROCESS.KINE.resamp(:,1))/100),...
%             length(PROCESS.KINE.resamp(:,1)));  % Head GCMA time vector  
     disp(['____ALL data resampled to ' num2str(fs) ' Hz'])
     
%% Precise alignment of both target vector and input data streams
    % Take out extra samples in the beginning. Should only be 1-5 samples to remove
     samplesrmvd = abs(length(PROCESS.EEG.resamp) - length(PROCESS.CLASS.classlabel));
    if samplesrmvd > 5;
        disp('Resampled data streams larger than five samples, re-check syncing alignment.');
    end
    [PROCESS.EEG.resamp,PROCESS.CLASS.classlabel] = aligntimesamplesize(PROCESS.EEG.resamp,PROCESS.CLASS.classlabel);

    disp('____target vector and EEG data time sample numbers aligned') 
    
%%  Filter EEG within a specific band
    n_f=3;  % filter order
    bpass_freq=[1 4]; % band-pass frequencies
    PROCESS.EEG.filtered = filter_data_bpass_NOCELL(PROCESS.EEG.resamp, fs, n_f, bpass_freq);    
%     PROCESS.KINE.filtered = filter_data_bpass(PROCESS.KINE.resamp, fs, n_f, delta_freq);
    disp('____EEG data filtered')
    
%% Minimal Redundancy - Maximum Relevance Channel Feature Selection
data_temp = PROCESS.EEG.filtered;   target_temp = PROCESS.CLASS.classlabel;
%{
num_featList = [100, 200, 300];
rankFeat = cell(1,length(num_featList));        mrmrMI = cell(1,length(num_featList));
for num_feat = 1:length(num_featList);
    [rankFeat{num_feat}, mrmrMI{num_feat}] = mrmr_mid_d( data_temp, target_temp, num_featList(num_feat));
end
%
topNfeatures = floor(0.7 * size(PROCESS.EEG.filtered,2));
[rankFeat, mrmrMI] = mrmr_mid_d( data_temp, target_temp, topNfeatures );

% add to data structure "PROCESS"
PROCESS.FEATSELECT.MRMRrankedFeat = rankFeat;
for rf = 1:length(rankFeat); PROCESS.FEATSELECT.MRMRrankedChnLabels{rf} = PROCESS.EEG.channelOrder_rmvdCHNS{rankFeat(rf)}; end
PROCESS.FEATSELECT.MRMR_MI = mrmrMI;
 %}   
%%  Lag-Based Feature Extraction
    % for a sampling frequency of 100 Hz, the lags will vary from 0 to -90 ms in steps of 10 ms
    lag = (0:9); max_lag=10; %maximum lag of -90 ms 
    PROCESS.EEG.featmat = timelag(PROCESS.EEG.filtered,lag);     % generate feature matrix
%     PROCESS.KINE.featmat = timelag(PROCESS.KINE.filtered,lag);     % generate feature matrix
    PROCESS.CLASS.classlabel=PROCESS.CLASS.classlabel(max_lag:end);     % truncate target vector to maintain same number of rows
    disp('____Feature matrix produced')
   
%%  Standardization across channels (Z-scores: subtracting mean from each data point and deviding by standard deviation )
    %tip: mean 'should be already' zero after filtering, check this if your filter is working correctly
    PROCESS.EEG.zscores_featmat = zscore(PROCESS.EEG.featmat);
%     PROCESS.KINE.zscores_featmat = zscore(PROCESS.KINE.featmat);
    disp('____Feature matrix standardized')
disp('--------------------------------------------------')
timecompute(toc)
disp('--------------------------------------------------')
end  %repeat for each infant
return
%% Assigning Labels to each Feature
    featLabels = struct(); % define data structure of feature labels
    for el = 1:length(lag);
        for chns = 1:num_EEGchns;
            featLabels( (chns*length(lag)) - (length(lag)-el) ).Chan = PROCESS.EEG.channelOrder_rmvdCHNS{chns};
            featLabels( (chns*length(lag)) - (length(lag)-el) ).Lag = lag(el);
        end
    end
    PROCESS.CLASS.featLabels = featLabels; % add to data structure "PROCESS"  

%% Minimal Redundancy - Maximum Relevance Feature Selection
%{
data_temp = PROCESS.EEG.zscores_featmat;   target_temp = PROCESS.CLASS.classlabel;
%{
num_featList = [100, 200, 300];
rankFeat = cell(1,length(num_featList));        mrmrMI = cell(1,length(num_featList));
for num_feat = 1:length(num_featList);
    [rankFeat{num_feat}, mrmrMI{num_feat}] = mrmr_mid_d( data_temp, target_temp, num_featList(num_feat));
end
%}
topNfeatures = floor(0.5 * size(PROCESS.EEG.zscores_featmat,2));
[rankFeat, mrmrMI] = mrmr_mid_d( data_temp, target_temp, topNfeatures );

% add to data structure "PROCESS"
PROCESS.FEATSELECT.MRMRrankedFeat = rankFeat;
for rf = 1:length(rankFeat); PROCESS.FEATSELECT.MRMRrankedFeatLabels{rf} = PROCESS.CLASS.featLabels(rankFeat(rf)); end
PROCESS.FEATSELECT.MRMR_MI = mrmrMI;
%}
   
%% LFDA-GMM classifier initializations

%---Classifier Input Definitions-----------------------------------------------%
FeatMatXsubj = PROCESS.EEG.zscores_featmat;  %data_temp( :,(rankFeat) );  % time samples-by-channels and their corresponding lags
[num_tsamples, num_features] = size(FeatMatXsubj); 
ClassLblXsubj = PROCESS.CLASS.classlabel; %target_temp;   % time samples-by-1 target vector 
%------------------------------------------------------------------------------%
%---Optimized LFDA List Initializations (Only when analyzing multiple subjects)%
% Infants ==>  N09 | J20 | B06 | GR09 | A06 | LW10 | AR16 | RB23
optdimList = [ 225,  225,  200,   225,  250,   175,   125,   225 ];
optknnList = [   7,    7,    7,    19,   19,    13,     9,    19 ];
%------------------------------------------------------------------------------%
%---Optimization Initializations-----------------------------------------------%
    dimIDX=1;  % for indexing the dimension parameter
    knnIDX=1;  % for indexing the kNN parameter
    prcntTrnTst = (2/3)*100; % train-test percentage value to use (validating is 1 - tr_perc)
    tpIDX=1;   % for indexing the percentage of training samples
    prcntTrain = (0.5)*100; % training percentage value to use (testing is 1 - tr_perc)   
    dimSet = (floor(0.25*num_features):floor(0.75*num_features)); %optdimList(infant); %  % set of values for optimization of LDFA parameter 'dim'
    knnSet = (7:2:21); %optknnList(infant);  %  % set of values for optimization of LDFA parameter 'knn'
    maxiterations = 10; % maximum number of iterations to run for the LFDA GMM classifier due to random sampling. More than 10 iterations preferred 
%------------------------------------------------------------------------------%


%% Re-arrangement of data by targeted classes
% Providing a percentage of training and testing samples to use for the classifier
    [data1,C1] = ytoc(FeatMatXsubj',ClassLblXsubj);
    
% Remove classes with very small sample sizes
    [ data1, C1 ] = removesmallNclasses( data1, C1 );
    
% define the numbers of behaviors (classes) used for classification
if C1; clear num_classes; num_classes = length(C1)-1; end
    
% save class labels 'classesSelected' into data structure 'PROCESS'
    PROCESS.CLASS.classesSelected = classesSelected;
    
% percentage of training-testing data samples to use (using percentage of least populated class as reference)    
    num_trntst=round(min(C1)*prcntTrnTst/100);
    
% number of validating data samples to use (smallest class size - number of training samples)
    num_valid=min(C1)-num_trntst;
    
% randomly split dataset into a set for training and testing the
% classifier as well as an unseen set to validate the classifier
    [data_trntst, data_valid, c_trntst, c_valid,...
        rnd_indices_trntst, rnd_indices_valid] = ...
        datasplitfcn(data1((C1(1):end),:), C1(2:end), '2', num_trntst, num_valid);
    
% percentage of training data samples to use (using percentage of any equivalent class as reference)    
    num_train=round(unique(c_trntst) * prcntTrain/100); 

% number of testing data samples to use (class size - number of training samples)
    num_test=unique(c_trntst) - num_train;
    
% perform LFDA-GMM classification through a grid search
[ ClassifierResults ] = LFDAGMMGridSearch_wParComp( data_trntst, c_trntst, num_train, num_test, dimSet, knnSet, maxiterations  );
    
%{    
% open Parallel Computing Toolbox
if matlabpool('SIZE') == 0;   matlabpool('open'); end
classmatrixlist_LFDAgridsearch = cell(length(dimSet),length(knnSet));
confusionmatrixlist_LFDAgridsearch = cell(length(dimSet),length(knnSet));
accuracylist_LFDAgridsearch = cell(length(dimSet),length(knnSet));
precisionlist_LFDAgridsearch = cell(length(dimSet),length(knnSet));
sensitivitylist_LFDAgridsearch = cell(length(dimSet),length(knnSet));
GMmodel_LFDAgridsearch = cell(length(dimSet),length(knnSet));
TransformMat_LFDAgridsearch = cell(length(dimSet),length(knnSet));
mean_accuracy_LFDAgridsearch = zeros(length(dimSet),length(knnSet));
std_accuracy_LFDAgridsearch = zeros(length(dimSet),length(knnSet));
mean_precision_LFDAgridsearch = cell(length(dimSet),length(knnSet));
std_precision_LFDAgridsearch = cell(length(dimSet),length(knnSet));
mean_sensitivity_LFDAgridsearch = cell(length(dimSet),length(knnSet));
std_sensitivity_LFDAgridsearch = cell(length(dimSet),length(knnSet));
mean_F1score_LFDAgridsearch = cell(length(dimSet),length(knnSet));
SummedConfusionMatrix = cell(length(dimSet),length(knnSet));
tic;
parfor dim = 1:length(dimSet);  % for optimization of LDFA parameter 'dim'
    disp(['Computing LFDA Grid Search on Parallel Computing Worker ', num2str(dim)])
%     knnIDX=1;  % for indexing the kNN parameter
% slicing variables for Parallel Computing 'FOR' loop
    classmatrixlist_dimSlice = cell(1,length(knnSet));
    confusionmatrixlist_dimSlice = cell(1,length(knnSet));
    accuracylist_dimSlice = cell(1,length(knnSet));
    precisionlist_dimSlice = cell(1,length(knnSet));
    sensitivitylist_dimSlice = cell(1,length(knnSet));
    GMmodel_dimSlice = cell(1,length(knnSet));
    TransformMat_dimSlice = cell(1,length(knnSet));
    mean_accuracy_dimSlice = zeros(1,length(knnSet));
    std_accuracy_dimSlice = zeros(1,length(knnSet));
    mean_precision_dimSlice = cell(1,length(knnSet));
    std_precision_dimSlice = cell(1,length(knnSet));
    mean_sensitivity_dimSlice = cell(1,length(knnSet));
    std_sensitivity_dimSlice = cell(1,length(knnSet));
    mean_F1score_dimSlice = cell(1,length(knnSet));
    SummedConfusionMatrix_dimSlice = cell(1,length(knnSet));

    for knn = 1:length(knnSet); % for optimization of LDFA parameter 'knn'
        % slicing variables for Parallel Computing 'FOR' loop
        classmatrixlist = cell(1, maxiterations);
        confusionmatrixlist = cell(1, maxiterations);
        accuracylist = zeros(1, maxiterations);
        precisionlist = cell(1, maxiterations);
        sensitivitylist = cell(1, maxiterations);
        GMmodel = cell(1, maxiterations);
        TransformMat = cell(1, maxiterations);
        
        for iteration = 1:maxiterations;   % for computing the mean and standard deviation classification accuracy rates
        %% Random Sub-Sampling Cross Validation              
            % randomly split dataset into samples set for training and testing 
            [data_train, data_test, c_train, c_test, rnd_indices_train,...
                rnd_indices_test] = datasplitfcn(data_trntst,...
                c_trntst, '2', num_train, num_test);
            % shift randomized index values to reflect indices of 'data1' 
            indices_train = zeros(c_train(1,1), length(c_train)); % preallocate
            indices_test = zeros(c_test(1,1), length(c_test));    % preallocate
            for c = 1:length(c_train)
                indices_train(:,c) = rnd_indices_train(:,c) + sum(c_trntst(1:c));
                indices_test(:,c) = rnd_indices_test(:,c) + sum(c_trntst(1:c));
            end
        %% Local Fisher's Discriminant Analysis (LFDA)
%             dim=12; knn=1 ; %knn=odd number  
            % prepare input arguements for LFDA function
            X = data_train';            Y = ctoy(c_train);
            display('LFDA')
                                                         %    dim:        dimensionality of reduced space (default: d)
            [T,Z] = lfda(X,Y,dimSet(dim),'plain',knnSet(knn));           %    T  : d x r  transformation matrix (Z=T'*X)
                                                         %    Z  : r x n  matrix of dimensionality reduced samples 
                                                         %              metric: type of metric in the embedding space (default: 'weighted')
                                                         %              'weighted'        --- weighted eigenvectors 
                                                         %              'orthonormalized' --- orthonormalized
                                                         %              'plain'           --- raw eigenvectors               
            TransformMat{iteration} = T; % save transformation matrix
%                 r = 1; classes = [2,3];    class_names = {'Imitate','Point'};
%                 lfda_display( X, Y, T, Z, r, classes, class_names )

         %% Gaussian Mixture Modeling (GMM)
            gmm_data_train = Z;     % training data for GMM
            gmm_data_test  = data_test * T;     % testing data for GMM
            max_cluster = 10;  % Maximum number of clusters considered by GMM
            reg = 1e-5;        % Regularization parameter used in GMM
            display('GMM')
            try
            [class_matrix, density, posterior, comp_weight, obj, value] = gmmclassifier(gmm_data_train,  gmm_data_test,  c_train,  c_test,  max_cluster,  reg);
            catch err
                   disp('Error within the GMM classifier');
                   disp(err.getReport)
                continue 
            end
            % Save Classification Results into Data Structure
            classmatrixlist{iteration}=class_matrix;
            confusion_matrix = confusionmatrix_withSpecificities(class_matrix);      % compute confusion matrix with specificity per class included
            confusionmatrixlist{iteration} = confusion_matrix;
            accuracylist(iteration) = confusion_matrix(end,end-1);        % vector used to save accuracies
            precisionlist{iteration}=transpose(confusion_matrix(end,1:end-2));
            sensitivitylist{iteration}=confusion_matrix(1:end-1,end-1);
            GMmodel{iteration}=obj;


            display(['dimension value at ',num2str(dimSet(dim)),', knn value at ',num2str(knnSet(knn))]);
            display(['Iteration number ',num2str(iteration),', and accuracy: ',...
                num2str(accuracylist(iteration))]);
        end % repeat for each iteration of the LFDA-GMM classifier (for random sub-sampling X-validation)
       
        
        % add saved variables from 'iteration' loop to slicing variable for
        % each 'knn' iteration
        classmatrixlist_dimSlice{knn} = classmatrixlist;
        confusionmatrixlist_dimSlice{knn} = confusionmatrixlist;
        accuracylist_dimSlice{knn} = accuracylist;
        precisionlist_dimSlice{knn} = precisionlist;
        sensitivitylist_dimSlice{knn} = sensitivitylist;
        GMmodel_dimSlice{knn} = GMmodel;
        TransformMat_dimSlice{knn} = TransformMat;
        
        % additional calculations to perform and results to save as slicing
        % variables for each 'knn' iteration
        mean_accuracy_dimSlice(knn) = mean(accuracylist);
        std_accuracy_dimSlice(knn) = std(accuracylist);
        mean_precision_dimSlice{knn} = mean([precisionlist{:}],2);
        std_precision_dimSlice{knn} = std([precisionlist{:}],0,2);
        mean_sensitivity_dimSlice{knn} = mean([sensitivitylist{:}],2);
        std_sensitivity_dimSlice{knn} = std([sensitivitylist{:}],0,2);
        mean_F1score_dimSlice{knn} = ...
            2*(([mean_precision_dimSlice{knn}].*[mean_sensitivity_dimSlice{knn}])./...
            (([mean_precision_dimSlice{knn}] + [mean_sensitivity_dimSlice{knn}])));
        SummedConfusionMatrix_dimSlice{knn} = confusionmatrix_withSpecificities([classmatrixlist{:}]);       
    end         % repeat for each 'kNN' value

    classmatrixlist_LFDAgridsearch(dim,:) = classmatrixlist_dimSlice;
    confusionmatrixlist_LFDAgridsearch(dim,:) = confusionmatrixlist_dimSlice;
    accuracylist_LFDAgridsearch(dim,:) = accuracylist_dimSlice;
    precisionlist_LFDAgridsearch(dim,:) = precisionlist_dimSlice;
    sensitivitylist_LFDAgridsearch(dim,:) = sensitivitylist_dimSlice;
    GMmodel_LFDAgridsearch(dim,:) = GMmodel_dimSlice;
    TransformMat_LFDAgridsearch(dim,:) = TransformMat_dimSlice;
    mean_accuracy_LFDAgridsearch(dim,:) = mean_accuracy_dimSlice;
    std_accuracy_LFDAgridsearch(dim,:) = std_accuracy_dimSlice;
    mean_precision_LFDAgridsearch(dim,:) = mean_precision_dimSlice;
    std_precision_LFDAgridsearch(dim,:) = std_precision_dimSlice;
    mean_sensitivity_LFDAgridsearch(dim,:) = mean_sensitivity_dimSlice;
    std_sensitivity_LFDAgridsearch(dim,:) = std_sensitivity_dimSlice;
    mean_F1score_LFDAgridsearch(dim,:) = mean_F1score_dimSlice;
    SummedConfusionMatrix(dim,:) = SummedConfusionMatrix_dimSlice;
end    % repeat for each 'dim' value
disp('--------------------------------------------------')
timecompute(toc)
disp('--------------------------------------------------')
matlabpool('close') % close Parallel Computing Toolbox
%}
%% Saving all classifier results to data structure
PROCESS.CLASSIFIER = ClassifierResults;
%{
PROCESS.CLASSIFIER.classmatrixlist = classmatrixlist_LFDAgridsearch;
PROCESS.CLASSIFIER.confusionmatrixlist = confusionmatrixlist_LFDAgridsearch;
PROCESS.CLASSIFIER.accuracylist = accuracylist_LFDAgridsearch;       
PROCESS.CLASSIFIER.precisionlist = precisionlist_LFDAgridsearch;
PROCESS.CLASSIFIER.sensitivitylist = sensitivitylist_LFDAgridsearch;
PROCESS.CLASSIFIER.GMmodel = GMmodel_LFDAgridsearch;
PROCESS.CLASSIFIER.TransformMat = TransformMat_LFDAgridsearch;
PROCESS.CLASSIFIER.mean_accuracy = mean_accuracy_LFDAgridsearch;
PROCESS.CLASSIFIER.std_accuracy =std_accuracy_LFDAgridsearch;
PROCESS.CLASSIFIER.mean_precision = mean_precision_LFDAgridsearch;
PROCESS.CLASSIFIER.std_precision = std_precision_LFDAgridsearch;
PROCESS.CLASSIFIER.mean_sensitivity = mean_sensitivity_LFDAgridsearch;
PROCESS.CLASSIFIER.std_sensitivity = std_sensitivity_LFDAgridsearch;
PROCESS.CLASSIFIER.mean_F1score = mean_F1score_LFDAgridsearch;
PROCESS.CLASSIFIER.SummedConfusionMatrix = SummedConfusionMatrix;
PROCESS.LFDA_OPTIMIZATION.rSet = dimSet;
PROCESS.LFDA_OPTIMIZATION.kNNSet = knnSet;
%}

for r=1:length(dimSet);
    for k=1:length(knnSet);
        PROCESS.LFDA_OPTIMIZATION.acc_mean(r,k)= PROCESS.CLASSIFIER.mean_accuracy(r,k);
        PROCESS.LFDA_OPTIMIZATION.acc_std(r,k) = PROCESS.CLASSIFIER.std_accuracy(r,k);
    end
end

%}
cd([serverPath1 '/Classification Results'])
diary('ON');  % open diary file used for saving classification info
disp(['Completion of classification analysis with low freq. (' ...
    num2str(bpass_freq(1)) '-' num2str(bpass_freq(2)) ' Hz) EEG data',...
    '[using top ', num2str(topNfeatures), ' features based on  MR-MR feature selection code]'])

% find the indices of the LFDA parameters that yield the highest mean accuracy
[max_row, max_col]=find(PROCESS.LFDA_OPTIMIZATION.acc_mean == max(PROCESS.LFDA_OPTIMIZATION.acc_mean(:)));

% ...and find the mean and standard deviation accuracy
opt_OAmn = PROCESS.LFDA_OPTIMIZATION.acc_mean(max_row(1),max_col(1));
opt_OAstd = PROCESS.LFDA_OPTIMIZATION.acc_std(max_row(1),max_col(1));

% ...as well as the LFDA parameters ('r' and 'knn')
opt_r = PROCESS.LFDA_OPTIMIZATION.rSet(max_row(1));
opt_knn = PROCESS.LFDA_OPTIMIZATION.kNNSet(max_col(1));

% display the results
disp(['Optimized LFDA parameters: ' 10 '.......r = ' num2str(opt_r) 10 ...
    '.......kNN = ' num2str(opt_knn) 10 '.......Overall Accuracy = ' ...
    num2str(opt_OAmn) ' (+/-) ' num2str(opt_OAstd)])
diary('OFF');  % close diary file used for saving classification info
disp('--------------------------------------------------')
timecompute(toc)
disp('--------------------------------------------------')

[ ClassifierValidationResults ] = LFDAGMMValidateGridSearch( data_valid, c_valid, optGMmodel, optTransformMat, maxiterations );
%{
%% Validation of Optimal LFDA parameters on unseen data
% LFDA-GMM (perform transformation of data 
% and testing on all 10 model iterations
% and average all accuracies)
GMmodel = cell(1,length(c_valid));  % preallocate
for iter = 1:maxiterations
    display(['iteration #', num2str(iter)])
    display('....Transforming validation data into LFDA sub-space')
    % LFDA
    % prepare input arguments for LFDA function
    % validating data for GMM                                         
    gmm_data_valid  = data_valid * PROCESS.CLASSIFIER.TransformMat{max_row,max_col}{iter};         
    
    GMmodel = PROCESS.CLASSIFIER.GMmodel{max_row,max_col}{iter};
    display('....Gaussian mixture model validating')
    [ classmat_valid, ~, ~ ] = gmmclassifier_test( GMmodel, gmm_data_valid, c_valid );
    % Save Classification Results into Data Structure
    PROCESS.VALIDATION.classmatrixlist(:,:,iter)=classmat_valid;
    confusion_matrix = confusionmatrix_withSpecificities(classmat_valid);      % compute confusion matrix with specificity per class included
    PROCESS.VALIDATION.confusionmatrixlist(:,:,iter) = confusion_matrix;
    PROCESS.VALIDATION.accuracylist(iter) = confusion_matrix(end,end-1);        % vector used to save accuracies
    PROCESS.VALIDATION.precisionlist(iter,:)=confusion_matrix(end,1:end-2);
    PROCESS.VALIDATION.sensitivitylist(:,iter)=confusion_matrix(1:end-1,end-1); 
end
% Save More Classification Results into Data Structure
PROCESS.VALIDATION.mean_accuracy = mean(PROCESS.VALIDATION.accuracylist);
PROCESS.VALIDATION.std_accuracy = std(PROCESS.VALIDATION.accuracylist);
PROCESS.VALIDATION.mean_precision = mean(PROCESS.VALIDATION.precisionlist);
PROCESS.VALIDATION.std_precision = std(PROCESS.VALIDATION.precisionlist,0,1);
PROCESS.VALIDATION.mean_sensitivity = transpose(mean(PROCESS.VALIDATION.sensitivitylist,2));
PROCESS.VALIDATION.std_sensitivity = transpose(std(PROCESS.VALIDATION.sensitivitylist,0,2));
PROCESS.VALIDATION.mean_F1score = 2*((PROCESS.VALIDATION.mean_precision.*PROCESS.VALIDATION.mean_sensitivity)./...
    ((PROCESS.VALIDATION.mean_precision + PROCESS.VALIDATION.mean_sensitivity)));
PROCESS.VALIDATION.SummedConfusionMatrix = confusionmatrix_withSpecificities(PROCESS.VALIDATION.classmatrixlist);
%}
diary('ON'); % open diary file used for saving classification info
disp('Validation of optimal LFDA parameters on unseen data');
disp(['.......Overall Accuracy = ' num2str(PROCESS.VALIDATION.mean_accuracy) ...
    ' (+/-) ' num2str(PROCESS.VALIDATION.std_accuracy)]);
diary('OFF');  % close diary file used for saving classification info
disp('--------------------------------------------------')
timecompute(toc)
disp('--------------------------------------------------')

%% Display final confusion matrix in figure
% tasklabels_shift = circshift(capitalize(unique(TASKSEGMENTS.Task(ismember(TASKSEGMENTS.TaskLabel,classesSelected)))),-1);
tasklabels = capitalize(unique(PROCESS.CLASS.TASKSEGMENTS.Task(ismember(TASKSEGMENTS.TaskLabel,classesSelected))));
% to calculate percentage of data points belonging to each class over total number of data points
for cl = 1:num_classes-1; tasklabels{cl,2} = (C1(cl+1)/sum(C1(2:end)))*100; end
% [ NormCM ] = ConfusionMatrix_display( PROCESS.CLASSIFIER(max_row,max_col).SummedConfusionMatrix, tasklabels );
[ NormCM ] = ConfusionMatrix_display( PROCESS.VALIDATION.SummedConfusionMatrix, tasklabels );
fig_hndl = figure(1);       enhanceFigure(fig_hndl); % enhance figure for export
mkdir([serverPath1,'\Classification Results'],'Plots') % make new 'Plots' folder
cd('Plots')
filename = ['ConfusionMatrix_',num2str(num_classes),'-class_',...
    num2str(bpass_freq(1)),'to',num2str(bpass_freq(2)),'Hz_wValidationStep_wMRMRfeatselection']; % name of saved file
saveas(gcf, [filename, '.fig'], 'fig'); % save as .fig file
hgexport(gcf, [filename, '.tif'], hgexport('factorystyle'), 'Format', 'tiff'); % save as .tif file       
cd('..') % go back to main folder
%% Saving the resulting variables as a MAT file
% keep only new EEG data
subfieldslist = {'comments','nbchan','trials','pnts','srate','xmin','xmax',...
    'times','data','chanlocs','chaninfo','ref','event','eventdescription',...
    'reject','gain','ASRdata','eventTriggers','ExperStartInd','ExperEndInd',...
    'uVdata','channelOrder','uVdata_rmvdCHNS','channelOrder_rmvdCHNS',...
    'timesampleStart','timesampleEnd','uVdata_syncd','uVdata_rmvdCHNS_syncd'}; 
for fields = 1:length(subfieldslist); PROCESS.EEG = rmfield(PROCESS.EEG,subfieldslist{fields}); end

% save data structure to 'classification results folder'
fields_to_save = {'EEG','CLASS','CLASSIFIER','FEATSELECT','LFDA_OPTIMIZATION','VALIDATION'};
filename = ['LFDA_Optimization_',num2str(num_classes),'-class_',...
    num2str(bpass_freq(1)),'to',num2str(bpass_freq(2)),'Hz_wValidationStep_wMRMRfeatselection.mat'];
cd([serverPath1 '\Classification Results']) % save to Infantdata path where 'Data' is located

save(filename, '-struct', 'PROCESS', fields_to_save{:});
disp('____EEG data structure saved to lab server folder "Data"')
disp('--------------------------------------------------')
timecompute(toc)
disp('--------------------------------------------------')
close all
% end   % repeat for each infant