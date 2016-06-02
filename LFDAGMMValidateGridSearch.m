function [ ClassifierValidationResults ] = LFDAGMMValidateGridSearch( data_valid, c_valid, optGMmodel, optTransformMat, maxiterations )
%% Validation of Optimal LFDA parameters on unseen data
% LFDA-GMM (perform transformation of data 
% and testing on all 10 model iterations
% and average all accuracies)

% Pre-initializations
if ~exist('maxiterations','var')
    maxiterations = 10;
end

GMmodel = cell(1,length(c_valid));  % preallocate to open object variable
classmatrixlist = cell(1, maxiterations);
confusionmatrixlist = cell(1, maxiterations);
accuracylist = zeros(1, maxiterations);
precisionlist = cell(1, maxiterations);
sensitivitylist = cell(1, maxiterations);
for iter = 1:maxiterations
    display(['iteration #', num2str(iter)])
    display('....Transforming validation data into LFDA sub-space')
    % LFDA
    % prepare input arguments for LFDA function
    % validating data for GMM                                         
    gmm_data_valid  = data_valid * optTransformMat{iter};         
    
    GMmodel = optGMmodel{iter};
    display('....Gaussian mixture model validating')
    [ classmat, ~, ~ ] = gmmclassifier_test( GMmodel, gmm_data_valid, c_valid );
    
    % Saving each measured classification result
    classmatrixlist{iter}=classmat;
    confusion_matrix = confusionmatrix_withSpecificities(classmat);      % compute confusion matrix with specificity per class included
    confusionmatrixlist{iter} = confusion_matrix;
    accuracylist(iter) = confusion_matrix(end,end-1);        % vector used to save accuracies
    precisionlist{iter}=confusion_matrix(end,1:end-2);
    sensitivitylist{iter}=confusion_matrix(1:end-1,end-1); 
end

% Calculation of measured classification results
mean_accuracy = mean(accuracylist);
std_accuracy = std(accuracylist);
mean_precision = mean([precisionlist{:}],2);
std_precision = std([precisionlist{:}],0,2);
mean_sensitivity = mean([sensitivitylist{:}],2);
std_sensitivity = std([sensitivitylist{:}],0,2);
mean_F1score = 2*((mean_precision.*mean_sensitivity)./...
    ((mean_precision + mean_sensitivity)));
SummedConfusionMatrix = confusionmatrix_withSpecificities([classmatrixlist{:}]);       

% Save into a main data structure
ClassifierValidationResults = struct(...
    'classmatrixlist',classmatrixlist,...
    'confusionmatrixlist',confusionmatrixlist,...
    'accuracylist',accuracylist,...
    'precisionlist',precisionlist,...
    'sensitivitylist',sensitivitylist,...
    'mean_accuracy',mean_accuracy,...
    'std_accuracy',std_accuracy,...
    'mean_precision',mean_precision,...
    'std_precision',std_precision,...
    'mean_sensitivity',mean_sensitivity,...
    'std_sensitivity',std_sensitivity,...
    'mean_F1score',mean_F1score,...
    'SummedConfusionMatrix',SummedConfusionMatrix);
end

