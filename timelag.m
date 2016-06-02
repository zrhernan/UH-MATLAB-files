% New, correct version of timelag 
% He, 07-01-13 

function output = timelag(EEG, lag) 
% Add time lag in EEG data. 
% [lag] is a row vector containing the lags we want to add, for instance lag = [10:10:50,100:100:500];
% time 0 is not included in result if there is no 0 entry in [lag]
% If lag==0, no lag is added. output is the same as input
% Example:
% a = [
% 
%     0.0300    0.5341    0.9397    0.2442
%     0.5357    0.8854    0.3545    0.2955
%     0.0871    0.8990    0.4106    0.6802
%     0.8021    0.6259    0.9843    0.5278
%     0.9891    0.1379    0.9456    0.4116
%     0.0669    0.2178    0.6766    0.6026
%     0.9394    0.1821    0.9883    0.7505
%     0.0182    0.0418    0.7668    0.5835
%     0.6838    0.1069    0.3367    0.5518
%     0.7837    0.6164    0.6624    0.5836];
%     
% b = timelag(a,[0 1 2]);
% b =
% 
%     0.0871    0.5357    0.0300    0.8990    0.8854    0.5341    0.4106    0.3545    0.9397    0.6802    0.2955    0.2442
%     0.8021    0.0871    0.5357    0.6259    0.8990    0.8854    0.9843    0.4106    0.3545    0.5278    0.6802    0.2955
%     0.9891    0.8021    0.0871    0.1379    0.6259    0.8990    0.9456    0.9843    0.4106    0.4116    0.5278    0.6802
%     0.0669    0.9891    0.8021    0.2178    0.1379    0.6259    0.6766    0.9456    0.9843    0.6026    0.4116    0.5278
%     0.9394    0.0669    0.9891    0.1821    0.2178    0.1379    0.9883    0.6766    0.9456    0.7505    0.6026    0.4116
%     0.0182    0.9394    0.0669    0.0418    0.1821    0.2178    0.7668    0.9883    0.6766    0.5835    0.7505    0.6026
%     0.6838    0.0182    0.9394    0.1069    0.0418    0.1821    0.3367    0.7668    0.9883    0.5518    0.5835    0.7505
%     0.7837    0.6838    0.0182    0.6164    0.1069    0.0418    0.6624    0.3367    0.7668    0.5836    0.5518    0.5835

if lag == 0
    output = EEG;
    return
end

lag = sort(lag,'ascend');

len = length(EEG)-lag(end);
output = zeros(len,size(EEG,2)*length(lag));
temp = zeros(length(lag),size(EEG,2));
disp('Constructing Feature Matrix....');
% TL = waitbar(0,'Constructing Feature Matrix....');
for i = 1: len
    for j = 1:length(lag)
        temp(j,:) = EEG(i+lag(end)-lag(j),:);
    end
    output(i,:) = reshape(temp,1,[]);
    disp(['Progress: ', num2str((i/len)*100),'%'])
%     waitbar(i/len);
end
% close(TL) %close the waitbar
end % EOF