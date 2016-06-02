%% Annotation for Infant Data
InfantDataAnnotList = {'N09-12-04-2013','J20-08-19-2013',...
    'B06-10-30-2013','GR09-07-12-2014','A06-09-28-2013','LW10-06-19-2014'};

EEGheaderList = {'nata-12-04-2013','john-08-19-2013','byran-10-30-2013',...
    'InfantProject_000001_GR_2014_07_12','autumn-09-28-2013',...
    'InfantProject_LW_06-19-14'};

%% Selection of infants
infant = menu(['Which infant would you',10,' like to analyze?'],'N09',...
    'J20','B06','GR09','A06','LW10');
close

%% Path directory initializations
serverPath = ['\\172.27.216.40\Contreras-UH\Infantdata\Data\',...
    InfantDataAnnotList{infant},'\EEG'];
desktopPath = ['C:\Users\zrhernan\Infant_decoding_files\Data\',...
    InfantDataAnnotList{infant}];
%% Opening the VHDR file
EEG = pop_loadbv(serverPath,[EEGheaderList{infant},'.vhdr']);
EEG.gain = 0.1;      % add gain for changing amplitude order of magnitude

%% Opening Event Markers (triggers) and plotting them as a bar graph
triggers = zeros(1,length(EEG.event));
for trig = 1:length(EEG.event)
    triggers(trig) = EEG.event(trig).latency;
end
bar(triggers)
EEG.triggers = triggers;

%% Input the Experiment Start and Stop Times

EEG.ExperStartInd = input('Based on the bar plot please input the index corresponding to the start of the experiment');
EEG.ExperEndInd = input('Based on the bar plot please input the index corresponding to the end of the experiment');

%% Saving the resulting variables as a MAT file
cd(serverPath) % save to Infantdata path
save('EEGfiles.mat', '-struct', 'EEG', 'comments', 'nbchan', 'trials',...
    'pnts', 'srate', 'xmin', 'xmax', 'times', 'data', 'chanlocs',...
    'chaninfo', 'ref', 'event', 'eventdescription', 'reject', 'gain', 'triggers');


cd(desktopPath) % save to path on Zach's desktop
save('EEGfiles.mat', '-struct', 'EEG', 'comments', 'nbchan', 'trials',...
    'pnts', 'srate', 'xmin', 'xmax', 'times', 'data', 'chanlocs',...
    'chaninfo', 'ref', 'event', 'eventdescription', 'reject', 'gain', 'triggers');