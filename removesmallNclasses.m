function [ data1, C1, classesSelected ] = removesmallNclasses( data1, C1, ClassSize_thr )
%   Check if the least populated class is too small, 
%       a.k.a. if the number of class samples is less than the reduced
%       number of dimensions,
%   and remove them from both the vector of N classes (C1) and the array
%   of EEG data (data1).
% Created by: Zachery Hernandez, University of Houston, 2015

% Place default class size threshold to be 50% of the smallest class size
if nargin < 3
    ClassSize_thr = round(min(C1)*50/100);
end

% to keep track of which classes were kept for further classification    
classesSelected = (0:length(C1)-1);
    
% while a percentage of the class size is too small...    
while size(data1,2) > ClassSize_thr
    % remove class number from vector of kept classes
    classesSelected(C1==min(C1)) = [];

    % create a start and...
    classMinI = sum(C1(1:find(C1==min(C1))-1));

    %...end index for removing data corresponding to the class
    classMinF = sum(C1(1:find(C1==min(C1))))-1;

    % and remove that section
    data1(classMinI:classMinF,:) = [];

    % remove the class size from vector of class sizes
    C1(C1==min(C1)) = [];
end

end  %EOF

