function [new_data,C] = ytoc(data,y)
% This function finds all the indices at which a particular class occurs
% at that sample and uses that to find the number of samples each class
% contains.
% INPUT:
%   data = tempmrx'        640 x alot
%   ic <-- all elements    640 x alot
% OUTPUT:
%   new_data = Sorted Data (index from 1 to n_end time samples)  
%   C = 1 x num_label vector of indices with the number of samples per class

y = y(:);   
num_label = length(unique(y));      %numlabel=4
class_label = unique(y);            %unique(y)=[1 2 3 4]'

ord_ran = [];
for i = 1 : num_label
    class_ind = find(y==class_label(i));
    C(i) = length(class_ind);
    ord_ran = [ord_ran;class_ind];
end

new_data = data(:,ord_ran)';