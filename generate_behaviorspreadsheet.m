%% for calling helper functions
addpath('\\172.27.216.40\Contreras-UH\Infantdata\Infantdata\Source-Estimation\uh_fieldtrip')
addpath(genpath('\\172.27.216.40\Contreras-UH\Infantdata\Infantdata\Source-Estimation\uh_fieldtrip\demos'))
addpath(genpath('\\172.27.216.40\Contreras-UH\Infantdata\Infantdata\Source-Estimation\uh_fieldtrip\fileio'))
addpath(genpath('\\172.27.216.40\Contreras-UH\Infantdata\Infantdata\Source-Estimation\uh_fieldtrip\preprocess'))
addpath(genpath('\\172.27.216.40\Contreras-UH\Infantdata\Infantdata\Source-Estimation\uh_fieldtrip\time_freq'))
addpath(genpath('\\172.27.216.40\Contreras-UH\Infantdata\Infantdata\Source-Estimation\uh_fieldtrip\template'))
addpath(genpath('\\172.27.216.40\Contreras-UH\Infantdata\Infantdata\Source-Estimation\uh_fieldtrip\utils'))
addpath(genpath('\\172.27.216.40\Contreras-UH\Infantdata\Infantdata\Source-Estimation\uh_fieldtrip\visualization'))
addpath(genpath('C:\Users\zrhernan\Documents\MATLAB\fieldtrip'))
addpath('\\bmi-nas-01\Contreras-UH\Infantdata\Infantdata\code\Zachs_Infant_decoding_files')

%% Generate List of Infant Data Folders
InfantDir = '\\172.27.216.40\Contreras-UH\Infantdata\Infantdata\Data\';
disp_flag = 0;
[InfantDataAnnotList,InfantID] = defineInfantFolders(InfantDir,disp_flag);

%% Extract behavior information
MainListofBehaviors = [];
for ii=1:length(InfantID)
     disp(['Opening EEG recording of subject ',InfantID{ii}])
    infantfolder = InfantDataAnnotList{ii};
    
    % infantfolder = InfantDataAnnotList{53}; % OL21
    serverPath1 = [InfantDir,infantfolder,'\EEG\'];
    cd(serverPath1)
    
    
    %%
    behavior_filename = [InfantDir,infantfolder,'\Behavioral Segmentation\class start-stop times.txt'];
    classinfo = uh_importRawEvent_wHandInfo(behavior_filename);
    if isempty(classinfo), disp('List of behaviors does not exist. Moving to next subject.'), continue, end
    %%
    disp('Appending behavioral list to main list...')
    subjIDlist = ['Subject ID'; repmat({InfantID{ii}},[size(classinfo,1)-1,1])];
    agelist = ['Age'; repmat( cellstr(InfantID{ii}(regexp(InfantID{ii},'\d'))), [size(classinfo,1)-1,1] )];
    classinfo = [subjIDlist,agelist,classinfo];
    if ii==1
        MainListofBehaviors = [MainListofBehaviors; classinfo];
    else
        MainListofBehaviors = [MainListofBehaviors; classinfo(2:end,:)];
    end
    
    clear classinfo
    
end

%% To change to a dataset
FullBehaviorListDS = cell2dataset(MainListofBehaviors);
FullBehaviorListDS.Age = str2double(FullBehaviorListDS.Age);
FullBehaviorListDS.StartTimeSample = str2double(FullBehaviorListDS.StartTimeSample);
FullBehaviorListDS.EndTimeSample = str2double(FullBehaviorListDS.EndTimeSample);
FullBehaviorListDS.TaskLevel = str2double(FullBehaviorListDS.TaskLevel);
FullBehaviorListDS.Action = lower(FullBehaviorListDS.Action);

%% sort by age
FullBehaviorListDS_Agesort = sortrows(FullBehaviorListDS,{'Age','Action'});

