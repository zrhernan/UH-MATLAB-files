function y = ctoy(C)
%------------------------------------------------------------------------------------
% Usage:
%   y = makeclassindex(C)
%
% Input:
%   C: vector of indices giving the number of samples per class
%
% Output:
%   y: n dimensional vector of class labels
%------------------------------------------------------------------------------------

y = [];
for i = 1 : length(C)
    y = [y; i*ones(C(i),1)];
end


