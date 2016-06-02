function [class_matrix, posterior, comp_weight, obj, value, minBIC] = gmmclassifier(data_train,data_test,c_train,c_test,max_cluster,reg)
%------------------------------------------------------------------------------------
% Gaussian Mixture Model Classifier
%
% Usage:
%   [class_matrix density posterior] = gmmclassifier(data_train,data_test,c_train,c_test,options)
%
% Input:
%   data_train:  n1 x d matrix of training samples
%                n1 -- number of training samples
%                d  -- number of original features
%   data_test:   n2 x d matrix of training samples
%                n2 -- number of testing samples
%                d  -- number of original features
%   c_train:     vector of indices giving number of training samples per class
%   c_test:      vector of indices giving number of training samples per class
%   max_cluster: maximum number of clusters to be considered
%   reg:         regularization parameter used in gaussian mixture model
%
% Output:
%   class_matrix: c x k matrix,
%                 where each entry is a number to which the test samples belongs
%                 c -- number of classes
%                 k -- number of samples in the largest class
%   posterior:    c x n2 matrix,
%                 where (i,j)-th entry is a j-th sample's posterior probability of class i
%
% Dependencies:
%   gmdistribution.fit.m
%
% Copyright 2012 --- Hyperspectral Image Analysis Lab, University of Houston
% http://hyperspectral.ee.uh.edu
%------------------------------------------------------------------------------------

num_c = length(c_train);
obj = cell(1,num_c);
cls_cluster_num = zeros(1,num_c);

ctr = 0;
for i = 1 : num_c
    temp_data = data_train((ctr+1):(c_train(i)+ctr),:);
    ctr = c_train(i) + ctr;
    model = cell(1,max_cluster);
    BIC = zeros(1,max_cluster);
    for j = 1 : max_cluster
        IDX = clusterdata(temp_data, 'maxclust', j);
        options = statset('MaxIter', 1000);
        model{j} = gmdistribution.fit(temp_data,j,'Start',IDX,...
            'CovType','full','Regularize',reg,'Options',options);
        BIC(j) = model{j}.BIC;
    end  % repeat to 'max_cluster'
    [minBIC, Components_num] = min(BIC);
   
    obj{i} = model{Components_num};
    comp_weight=obj{i}.PComponents;
    cls_cluster_num(i) = Components_num;
end % repeat to 'num_c'

[num_test, ~] = size(data_test);
density = zeros(num_c,num_test);
% class_matrix = zeros(num_c, num_test);
% value = zeros(num_c, num_test);

% This version is Very Slow !!!
% for i = 1 : num_test
%     for j = 1 : num_c
%         density(i,j) = pdf(obj{j},data_test(i,:));
%     end
%     [value, class(i)] = max(density(i,:));
% end
% posterior = density;

cumsum_c = [0 cumsum(c_test)];
for i = 1 : length(c_test)
    ctr = 0;
    for j = cumsum_c(i)+1 : cumsum_c(i+1)
        ctr = ctr + 1;
        for z = 1 : num_c
            density(z,j) = pdf(obj{z},data_test(j,:));
        end
        [value(ctr,i), class_matrix(ctr, i)] = max(density(:,j));
        
    end
end

class_matrix = class_matrix';
posterior = density;


