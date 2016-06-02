function [data_train,data_test,c_train,c_test,indices_train,indices_test] = datasplitfcn(data,C,type,num_train,num_test)
%------------------------------------------------------------------------------------
% Data Split
% 
% Usage:
%   [data_train,c_train,data_test,c_test] = datasplitfcn(data,C,options)
%
% Input:
%   data:      n x d matrix of samples
%              n -- number of total samples
%              d -- number of original features
%   C:         vector of indices giving the number of samples per class
%   type:      type of method to split the data
%   num_train: number (percentage) of training samples to choose
%   num_test:  number of testing samples to choose
%
% Output:
%   data_train: n1 x d matrix of training samples
%               n1 -- number of training samples
%   data_test:  n2 x d matrix of testing samples
%               n2 -- number of testing samples
%   c_train:    vector of indices giving the number of training samples per class
%   c_test:     vector of indices giving the number of testing samples per class
%------------------------------------------------------------------------------------

if nargin<2
    error('Please provide enough arguments!')
end

if nargin<3
  type = '1';
end

if nargin<4
    error('Please provide number of training samples to choose!')
end

if nargin<5
  num_test = [];
end

cumsum_c = [0,cumsum(C)];
num_c = length(C);
[n d] = size(data);
switch type
    % Split dataset based on number of training samples per class
    case '1'
        c_train = num_train * ones(1,num_c);
        c_test = C - num_train;
        cumsum_c1 = [0,cumsum(c_train)];
        cumsum_c2 = [0,cumsum(c_test)];
        data_train = zeros(sum(c_train),d);
        data_test = zeros(sum(c_test),d);
        for i = 1 : num_c
            tempdata = data(cumsum_c(i)+1:cumsum_c(i+1),:);
            rand_num = randperm(C(i));
            data_train(cumsum_c1(i)+1:cumsum_c1(i+1),:) = tempdata(rand_num(1:c_train(i)),:);
            data_test(cumsum_c2(i)+1:cumsum_c2(i+1),:) = tempdata(rand_num(c_train(i)+1:end),:);
        end 
    % Split dataset based on number of training and testing samples per class
    case '2'
        data_train = [];
        data_test = [];
        for i = 1 : num_c
            tempdata = data(cumsum_c(i)+1:cumsum_c(i+1),:);
            rand_num = randperm(C(i));
            indices_train(:,i) = rand_num(1:num_train);
            indices_test(:,i) = rand_num(num_train+1:num_train+num_test);
            data_train = [data_train;tempdata(indices_train(:,i),:)];
            data_test = [data_test;tempdata(indices_test(:,i),:)];
        end
        c_train = num_train .* ones(1,num_c);
        c_test = num_test .* ones(1,num_c);
    % Split dataset based on percentage of training samples per class.
    case '3'
        percent_train = num_train;
        c_train = ceil(percent_train * C);
        c_test = C - c_train;
        cumsum_c1 = [0,cumsum(c_train)];
        cumsum_c2 = [0,cumsum(c_test)];
        data_train = zeros(sum(c_train),d);
        data_test = zeros(sum(c_test),d);
        for i = 1 : num_c
            tempdata = data(cumsum_c(i)+1:cumsum_c(i+1),:);
            rand_num = randperm(C(i));
            data_train(cumsum_c1(i)+1:cumsum_c1(i+1),:) = tempdata(rand_num(1:c_train(i)),:);
            data_test(cumsum_c2(i)+1:cumsum_c2(i+1),:) = tempdata(rand_num(c_train(i)+1:end),:);
        end
end

