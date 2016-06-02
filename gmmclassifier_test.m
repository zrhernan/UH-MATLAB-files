function [ class_matrix, posterior, value ] = gmmclassifier_test( obj, data_test, c_test )
%------------------------------------------------------------------------------------
% Gaussian Mixture Model Classifier
%
% Usage:
%   [class_matrix density posterior value] = gmmclassifier_test(obj data_test c_test)
%
% Input:
%   obj:         1 x c Gaussian mixture distribution cell array
%                c  -- number of classes
%   data_test:   n2 x d matrix of training samples
%                n2 -- number of testing samples
%                d  -- number of original features
%   c_test:      vector of indices giving number of training samples per class
%
% Output:
%   class_matrix: c x k matrix,
%                 where each entry is a number to which the test samples belongs
%                 c -- number of classes
%                 k -- number of samples in the largest class
%   posterior:    c x n2 matrix,
%                 where (i,j)-th entry is a j-th sample's posterior probability of class i
%   value:        c x k matrix,
%                 where each entry is the maximum posterior probability value of a given test sample
%                 c -- number of classes
%                 k -- number of samples in the largest class
%
% Zachery Hernandez, Non-invasive Brain-Machine Interface Systems Lab, 2015
% code taken from Hyperspectral Image Analysis Lab, University of Houston
% http://hyperspectral.ee.uh.edu
%------------------------------------------------------------------------------------
num_c = length(c_test);
[num_test, ~] = size(data_test);
density = zeros(num_c,num_test);
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

end

