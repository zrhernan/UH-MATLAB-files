function plotEEGraster( EEGdata, options)
%PLOTEEGRASTER Generates a raster plot of EEG channels.
%INPUT: 
%   EEGdata --> array of sample-by-channel EEG time-series data. Can also
%               contain multiple arrays for channel comparison (sample-by-
%               channel-by-EEG array).                           
%   options --> contains the following optional attributes for the plot:                      
%       options.fighndl    --> figure handle
%       options.fs         --> sampling frequency (default = 1000 Hz)
%       options.ts         --> time sample values
%       options.chnlabels  --> EEG channel labels (default = channel numbers)
%       options.chanoffset --> EEG channel labels (default = 50)
%       options.clrset --> color to plot for EEG signals (default = {'k',
%                          'r','b','g'})
% Created by: Zach Hernandez, University of Houston, 2016
if nargin < 1
    error('No EEG data found.')
end

disp( 'Plotting EEG...' );

if ~exist('options','var')
    options = struct('fighndl',figure(111),...
        'fs',1000,...
        'ts',(0:1/fs:size(EEGdata,1)),...
        'chnlabels',{cellstr(num2str((1:size(EEGdata,2))'))},...
        'clrset',{{'k'}},...
        'chanoffset',10);
end

if ~isfield(options, 'fighndl'), options.fighndl = figure(111); end
if ~isfield(options, 'fs'), options.fs = 1000; end
if ~isfield(options, 'ts'), options.ts = (0:1/fs:size(EEGdata,1)); end
if ~isfield(options, 'chnlabels'), 
    options.chnlabels = cellstr(num2str((1:size(EEGdata,2))'));
end
if ~isfield(options, 'chanoffset'), options.chanoffset = 50; end
if ~isfield(options, 'clrset'), options.clrset = {'k','r','b','g'}; end

if length(EEGdata) ~= size(EEGdata,1), EEGdata = permute(EEGdata,[2 1 3]); end

[N,numchns,numsets] = size(EEGdata);
set(options.fighndl,'PaperPosition',[0.25 2.5 8 6]);
chnOS = options.chanoffset;

for dim3 = 1:numsets %for each EEG set to overlay (if more than 1)
    for chNo = 1:numchns,
        plot( EEGdata(:,chNo,dim3)-((chNo-1)*chnOS) , options.clrset{dim3});
        hold on;
    end
    lastYtick = -(numchns-1)*chnOS;
    set(gca,'xtick',(1:options.fs:N+1),'xticklabel',...
        options.ts(1:options.fs:N),'xlim',[0 N],...
        'ytick',(lastYtick:chnOS:0),'yticklabel',...
        flipud(options.chnlabels),'ylim',[lastYtick-chnOS chnOS],...
        'fontweight','bold','fontsize',12);
    xlabel('Time (s)','fontweight','bold','fontsize',20);
    ylabel('Channels','fontweight','bold','fontsize',20);
    hold off;
end

end    %EOF

