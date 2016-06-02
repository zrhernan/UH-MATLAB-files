function [featmat_new,classlabel_new] = aligntimesamplesize(featmat_old,classlabel_old)
%Aligns the time sample size of two time series data by truncating the
%longer time series
% WARNING: Should only be used for data with a few (< 5) samples difference

if length(classlabel_old) > size(featmat_old,1)
    sizediff = length(classlabel_old) - size(featmat_old,1);
    classlabel_new = classlabel_old(1:end-sizediff);
    featmat_new = featmat_old;
elseif size(featmat_old,1) > length(classlabel_old)
    sizediff = size(featmat_old,1) - length(classlabel_old);
    featmat_new = featmat_old(1:end-sizediff,:);
    classlabel_new = classlabel_old;
elseif size(featmat_old,1) == length(classlabel_old)
    featmat_new = featmat_old;
    classlabel_new = classlabel_old;
end

end %EOF