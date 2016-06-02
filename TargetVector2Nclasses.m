function [new_data,Nclasses] = TargetVector2Nclasses(data,TargetVector)
% This function finds all the indices at which a particular class occurs
% at that sample and uses that to find the number of samples each class
% contains.
% INPUT:
%          data : attributes x instances
%   TargetVector: 1 x instances
% OUTPUT:
%       new_data: Sorted Data (index from 1 to n_end time samples)  
%       Nclasses: 1 x num_label vector of indices with the number of samples per class

TargetVector = TargetVector(:);   
num_label = length(unique(TargetVector));      %numlabel=4
class_label = unique(TargetVector);            %unique(y)=[1 2 3 4]'

ord_ran = [];
for i = 1 : num_label
    class_ind = find(TargetVector==class_label(i));
    Nclasses(i) = length(class_ind);
    ord_ran = [ord_ran;class_ind];
end

new_data = data(:,ord_ran)';