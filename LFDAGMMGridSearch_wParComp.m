function [ ClassifierGridSearchResults ] = LFDAGMMGridSearch_wParComp( data_trntst, c_trntst, num_train, num_test, dimSet, knnSet, maxiterations  )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% Pre-initializations
if ~exist('dimSet','var')
    [~, num_features] = size(data_trntst);
    dimSet = (floor(0.25*num_features):300:floor(0.75*num_features));
end
if ~exist('knnSet','var')
    knnSet = (7:2:21);
end
if ~exist('maxiterations','var')
    maxiterations = 10;
end

% open Parallel Computing Toolbox
if isempty(gcp('nocreate')) == 1;   parpool('local'); end

% variables used in parallel computation
classmatrixlist = cell(length(dimSet),length(knnSet));
confusionmatrixlist = cell(length(dimSet),length(knnSet));
accuracylist = cell(length(dimSet),length(knnSet));
precisionlist = cell(length(dimSet),length(knnSet));
sensitivitylist = cell(length(dimSet),length(knnSet));
GMmodel = cell(length(dimSet),length(knnSet));
GMcompweights = cell(length(dimSet),length(knnSet));
GMposterior = cell(length(dimSet),length(knnSet));
classmatrixvaluelist = cell(length(dimSet),length(knnSet));
TransformMat = cell(length(dimSet),length(knnSet));
mean_accuracy = zeros(length(dimSet),length(knnSet));
std_accuracy = zeros(length(dimSet),length(knnSet));
mean_precision = cell(length(dimSet),length(knnSet));
std_precision = cell(length(dimSet),length(knnSet));
mean_sensitivity = cell(length(dimSet),length(knnSet));
std_sensitivity = cell(length(dimSet),length(knnSet));
mean_F1score = cell(length(dimSet),length(knnSet));
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
    GMcompweights_dimSlice = cell(1,length(knnSet));
    GMposterior_dimSlice = cell(1,length(knnSet));
    classmatrixvaluelist_dimSlice = cell(1,length(knnSet));
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
        classmatrixlist_knnSlice = cell(1, maxiterations);
        confusionmatrixlist_knnSlice = cell(1, maxiterations);
        accuracylist_knnSlice = zeros(1, maxiterations);
        precisionlist_knnSlice = cell(1, maxiterations);
        sensitivitylist_knnSlice = cell(1, maxiterations);
        GMmodel_knnSlice = cell(1, maxiterations);
        GMcompweights_knnSlice = cell(1, maxiterations);
        GMposterior_knnSlice = cell(1, maxiterations);
        classmatrixvaluelist_knnSlice = cell(1, maxiterations);
        TransformMat_knnSlice = cell(1, maxiterations);
        
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
            X = data_train';            Y = Nclasses2TargetVector(c_train);
            display('LFDA')
                                                         %    dim:        dimensionality of reduced space (default: d)
            [T,Z] = lfda(X,Y,dimSet(dim),'plain',knnSet(knn));           %    T  : d x r  transformation matrix (Z=T'*X)
                                                         %    Z  : r x n  matrix of dimensionality reduced samples 
                                                         %              metric: type of metric in the embedding space (default: 'weighted')
                                                         %              'weighted'        --- weighted eigenvectors 
                                                         %              'orthonormalized' --- orthonormalized
                                                         %              'plain'           --- raw eigenvectors               
            TransformMat_knnSlice{iteration} = T; % save transformation matrix
%                 r = 1; classes = [2,3];    class_names = {'Imitate','Point'};
%                 lfda_display( X, Y, T, Z, r, classes, class_names )

         %% Gaussian Mixture Modeling (GMM)
            gmm_data_train = Z;     % training data for GMM
            gmm_data_test  = data_test * T;     % testing data for GMM
            max_cluster = 10;  % Maximum number of clusters considered by GMM
            reg = 1e-5;        % Regularization parameter used in GMM
            display('GMM')
            try
            [class_matrix, posterior, comp_weight, obj, value] = gmmclassifier(gmm_data_train,  gmm_data_test,  c_train,  c_test,  max_cluster,  reg);
            catch err
                   disp('Error within the GMM classifier');
                   disp(err.getReport)
                continue 
            end
            % Save Classification Results into Data Structure
            classmatrixlist_knnSlice{iteration}=class_matrix;
            confusion_matrix = confusionmatrix_withSpecificities(class_matrix);      % compute confusion matrix with specificity per class included
            confusionmatrixlist_knnSlice{iteration} = confusion_matrix;
            accuracylist_knnSlice(iteration) = confusion_matrix(end,end-1);        % vector used to save accuracies
            precisionlist_knnSlice{iteration}=transpose(confusion_matrix(end,1:end-2));
            sensitivitylist_knnSlice{iteration}=confusion_matrix(1:end-1,end-1);
            GMmodel_knnSlice{iteration}=obj;
            GMcompweights_knnSlice{iteration}=comp_weight;
            GMposterior_knnSlice{iteration}=posterior;
            classmatrixvaluelist_knnSlice{iteration}=value;


            display(['dimension value at ',num2str(dimSet(dim)),', knn value at ',num2str(knnSet(knn))]);
            display(['Iteration number ',num2str(iteration),', and accuracy: ',...
                num2str(accuracylist_knnSlice(iteration))]);
        end % repeat for each iteration of the LFDA-GMM classifier (for random sub-sampling X-validation)
       
        
        % add saved variables from 'iteration' loop to slicing variable for
        % each 'knn' iteration
        classmatrixlist_dimSlice{knn} = classmatrixlist_knnSlice;
        confusionmatrixlist_dimSlice{knn} = confusionmatrixlist_knnSlice;
        accuracylist_dimSlice{knn} = accuracylist_knnSlice;
        precisionlist_dimSlice{knn} = precisionlist_knnSlice;
        sensitivitylist_dimSlice{knn} = sensitivitylist_knnSlice;
        GMmodel_dimSlice{knn} = GMmodel_knnSlice;
        GMcompweights_dimSlice{knn} = GMcompweights_knnSlice;
        GMposterior_dimSlice{knn} = GMposterior_knnSlice;
        classmatrixvaluelist_dimSlice{knn} = classmatrixvaluelist_knnSlice;
        TransformMat_dimSlice{knn} = TransformMat_knnSlice;
        
        % additional calculations to perform and results to save as slicing
        % variables for each 'knn' iteration
        mean_accuracy_dimSlice(knn) = mean(accuracylist_knnSlice);
        std_accuracy_dimSlice(knn) = std(accuracylist_knnSlice);
        mean_precision_dimSlice{knn} = mean([precisionlist_knnSlice{:}],2);
        std_precision_dimSlice{knn} = std([precisionlist_knnSlice{:}],0,2);
        mean_sensitivity_dimSlice{knn} = mean([sensitivitylist_knnSlice{:}],2);
        std_sensitivity_dimSlice{knn} = std([sensitivitylist_knnSlice{:}],0,2);
        mean_F1score_dimSlice{knn} = ...
            2*(([mean_precision_dimSlice{knn}].*[mean_sensitivity_dimSlice{knn}])./...
            (([mean_precision_dimSlice{knn}] + [mean_sensitivity_dimSlice{knn}])));
        SummedConfusionMatrix_dimSlice{knn} = confusionmatrix_withSpecificities([classmatrixlist_knnSlice{:}]);       
    end         % repeat for each 'kNN' value

    classmatrixlist(dim,:) = classmatrixlist_dimSlice;
    confusionmatrixlist(dim,:) = confusionmatrixlist_dimSlice;
    accuracylist(dim,:) = accuracylist_dimSlice;
    precisionlist(dim,:) = precisionlist_dimSlice;
    sensitivitylist(dim,:) = sensitivitylist_dimSlice;
    GMmodel(dim,:) = GMmodel_dimSlice;
    GMcompweights(dim,:) = GMcompweights_dimSlice;
    GMposterior(dim,:) = GMposterior_dimSlice;
    classmatrixvaluelist(dim,:) = classmatrixvaluelist_dimSlice;
    TransformMat(dim,:) = TransformMat_dimSlice;
    mean_accuracy(dim,:) = mean_accuracy_dimSlice;
    std_accuracy(dim,:) = std_accuracy_dimSlice;
    mean_precision(dim,:) = mean_precision_dimSlice;
    std_precision(dim,:) = std_precision_dimSlice;
    mean_sensitivity(dim,:) = mean_sensitivity_dimSlice;
    std_sensitivity(dim,:) = std_sensitivity_dimSlice;
    mean_F1score(dim,:) = mean_F1score_dimSlice;
    SummedConfusionMatrix(dim,:) = SummedConfusionMatrix_dimSlice;
end    % repeat for each 'dim' value
disp('--------------------------------------------------')
timecompute(toc)
disp('--------------------------------------------------')
delete(gcp) % close Parallel Computing Toolbox

% place all variables into a data structure
ClassifierGridSearchResults = struct(...
    'classmatrixlist',classmatrixlist,...
    'confusionmatrixlist',confusionmatrixlist,...
    'accuracylist',accuracylist,...
    'precisionlist',precisionlist,...
    'sensitivitylist',sensitivitylist,...
    'GMmodel',GMmodel,...
    'GMcompweights',GMcompweights,...
    'GMposterior',GMposterior,...
    'classmatrixvaluelist',classmatrixvaluelist,...
    'TransformMat',TransformMat,...
    'mean_accuracy',mean_accuracy,...
    'std_accuracy',std_accuracy,...
    'mean_precision',mean_precision,...
    'std_precision',std_precision,...
    'mean_sensitivity',mean_sensitivity,...
    'std_sensitivity',std_sensitivity,...
    'mean_F1score',mean_F1score,...
    'SummedConfusionMatrix',SummedConfusionMatrix);

end  %EOF

