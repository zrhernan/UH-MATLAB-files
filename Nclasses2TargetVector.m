function TargetVector = Nclasses2TargetVector(Nclasses)
%------------------------------------------------------------------------------------
% Usage:
%   TargetVector = Nclasses2TargetVector(Nclasses)
%
% Input:
%   Nclasses: vector of samples sizes per class
%
% Output:
%   TargetVector: n dimensional vector of class labels
%------------------------------------------------------------------------------------

TargetVector = [];
for i = 1 : length(Nclasses)
    TargetVector = [TargetVector; i*ones(Nclasses(i),1)];
end


