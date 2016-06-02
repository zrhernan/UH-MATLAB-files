function [ ordered_channel_list_chans_rmvd, varargout ] = removeEEGchannels( ordered_channel_list, chans_to_remove, varargin )
%REMOVEEEGCHANNELS Summary of this function goes here
%   Detailed explanation goes here

%% Optional Input Arguments
OIAflag = 0; % initialize optional input argument (OIA) flag
varargout{1} = [];
if exist('varargin','var')
    if strcmp(varargin{1},'EEG data');
        EEGdata = varargin{2};
        OIAflag = 1;
    else
        disp('Error:  Last input argument should be in the form of ''EEG data'', [name of EEG data variable]')
        return
    end
end

%% Use the 'chans_to_remove' to remove channels from 'channel_list'
chnsremoved=zeros(length(chans_to_remove),1);
for rmv = 1:length(chans_to_remove); 
    chnsremoved(rmv)=find(strcmp(ordered_channel_list, chans_to_remove{rmv})==1); 
end

%% Removing channels from the ordered list of channels
ordered_channel_list_chans_rmvd = ordered_channel_list;  % initialize channel order list of all channels
ordered_channel_list_chans_rmvd(:,chnsremoved)=[];   % and remove undesirable channels from the EEG channels list

%% Optional: Removing channels from the EEG data
% nargoutchk(0,1)
if OIAflag == 1;
    %first check if number of EEG data channels matches the 'ordered_channel_list'
    if size(EEGdata,2) ~= length(ordered_channel_list);
        disp('Error: ''EEGdata'' should contain the same number of channels as ''ordered_channel_list''')
        return
    end

    % now remove the channel from the EEG data array
    EEGdata_chans_rmvd = EEGdata;        % initialize EEG dataset of all channels
    EEGdata_chans_rmvd(:,chnsremoved)=[];   % and remove undesirable channels from the EEG dataset
    varargout{1} = EEGdata_chans_rmvd;
else
    return
end

end

