% clc, clear all, close all
tic;  % start computation time

% for calling helper functions
addpath('\\172.27.216.40\Contreras-UH\Infantdata\Infantdata\code\Zachs_Infant_decoding_files') 

%% Initialization of Data Structures
TimeLengthData = struct('durations',{},'mean',{},'std',{},'max',{},'min',{});

%% Annotation for Infant Data
InfantDataAnnotList = {'N09-12-04-2013','J20-08-19-2013','B06-10-30-2013',...
    'GR09-07-12-2014','A06-09-28-2013','LW10-06-19-2014','AR16-07-14-2014',...
    'RB23-12-04-2014','A18-08-15-2013'};

%% Selection of infants
InfantID = {'N09','J20','B06','GR09','A06','LW10','AR16','RB23','A18'};
month = [9, 20, 6, 9, 6, 10, 16, 23, 18];
TimeLengthData(1).BehaviorName = {'explore','imitate','observe','reach-grasp','reach-offer','rest'};
labelsTLengthPerInfant = {};        monthsPerBehavior = [];
meanTLengthPerInfant = [];          stdTLengthPerInfant = [];
for infant = 1:length(InfantID);
    %% Initializing and Assigning directory paths....
    serverPath = ['\\bmi-nas-01\Contreras-UH\Infantdata\Infantdata\code\Zachs_Infant_decoding_files\Data\',InfantDataAnnotList{infant}];
    %desktopPath = ['C:\Users\zrhernan\Infant_decoding_files\Data\',InfantDataAnnotList{infant}];
%     mkdir(serverPath,'Behavioral Segmentation')  
    cd(serverPath)

    %% Importation of List of Time-Segmented Tasks
    LastColumnIndex = [104,144,68,114,49,50,153,238,83];       LastColumnIndex_2class = [39,46,19]; %for binary class dataset
    LastColumnIndex_7C_J20 = 162; 
    TASKSEGMENTS = importTaskTrialInfo(['Data\',InfantDataAnnotList{infant},'\class start-stop times.txt'],...
        2, LastColumnIndex(infant)); %LastColumnIndex_7C_J20); %     % load list of start and stop times for each class
    % TASKSEGMENTS = importTaskTrialInfo(['Data\',InfantDataAnnotList{infant},'\Observe&Imitate_start-stop times.txt'],...
    %     2, LastColumnIndex_2class);  % load list of start and stop times for each 'observe' and 'imitate' class
%     PROCESS.CLASS.TASKSEGMENTS = TASKSEGMENTS;   % and save into PROCESS structure
    disp(['____List of tasks imported for ', InfantID{infant}])
    
    %% Calculate the time length + mean + standard deviation of each behavior
    BehaviorLabel = TimeLengthData.BehaviorName;
    for cl = 1:length(BehaviorLabel);
%         TimeLengthData.durations{1,cl} = unique(TASKSEGMENTS.Task(TASKSEGMENTS.TaskLabel==cl));
        TimeLengthData.durations{infant,cl} = TASKSEGMENTS.StopTime(strcmp(TASKSEGMENTS.Task,BehaviorLabel{cl}))...
            - TASKSEGMENTS.StartTime(strcmp(TASKSEGMENTS.Task,BehaviorLabel{cl}));
        TimeLengthData.mean(infant,cl) = mean(TimeLengthData.durations{infant,cl});
%         TimeLengthData.meanperinfant(infant,cl) = mean(TimeLengthData.durations{infant+1,:});
        if ~isempty(TimeLengthData.durations{infant,cl})
                TimeLengthData.std(infant,cl) = std(TimeLengthData.durations{infant,cl});
                TimeLengthData.max(infant,cl) = nanmax(TimeLengthData.durations{infant,cl});
                TimeLengthData.min(infant,cl) = nanmin(TimeLengthData.durations{infant,cl});
                TimeLengthData.InfantID{infant} = InfantID{infant};   
        end
    end
    TimeLengthData.durationsByInfant{infant} = vertcat(TimeLengthData.durations{infant,:});
    TimeLengthData.meanDurationsbyInfant(infant) = mean(TimeLengthData.durationsByInfant{infant});
    TimeLengthData.stdDurationsbyInfant(infant) = std(TimeLengthData.durationsByInfant{infant});
    %{
    % save to Infantdata path where 'Data' is located
    cd([serverPath '/Behavioral Segmentation'])
    % create filename
    filename = 'TimeLengthsandRelatedStatistics.mat';
    save(filename, '-struct', 'TimeLengthData');
    disp('____TimeLength data structure saved to lab server folder "Data"')

    %}
    
%     monthsPerBehavior = [monthsPerBehavior, repmat(month(infant),1,6)];
%     labelsTLengthPerInfant = {labelsTLengthPerInfant, TimeLengthData.durations{1,:}};
%     meanTLengthPerInfant = [meanTLengthPerInfant, TimeLengthData.mean];
%     stdTLengthPerInfant = [stdTLengthPerInfant, TimeLengthData.std];
end
disp('--------------------------------------------------')
timecompute(toc)
disp('--------------------------------------------------')


%% Plot the results
% gscatter(monthsPerBehavior,meanTLengthPerInfant,numericLabel)
errorbar(month, TimeLengthData.meanDurationsbyInfant/1000,...
    TimeLengthData.stdDurationsbyInfant/1000,'.r',...
    'markersize', 24, 'linewidth', 2)
ylabel('Average Behavior Duration (seconds)',...
            'FontSize',14,'FontWeight','bold')
xlabel('Age (months)','FontSize',14,'FontWeight','bold')
set(gca,'FontSize',14,'FontWeight','bold')
return


fig_hndl = figure;
for cl = 1:max(numericLabel);
    subplot(2,3,cl)
%     scatter(monthsPerBehavior(numericLabel==cl),...
%         (meanTLengthPerInfant(numericLabel==cl)/6000),80,'r','Filled')
    errorbar(monthsPerBehavior(numericLabel==cl),...
        (meanTLengthPerInfant(numericLabel==cl)/60000),...
        (stdTLengthPerInfant(numericLabel==cl)/60000),'.r',...
        'markersize',24,'linewidth',2)
    title(capitalize(TimeLengthData(1).durations{1,cl}),'FontSize',16,'FontWeight','bold')
    if ~mod((cl-1),3); 
        ylabel('Mean Time Length of Each Behavior (minutes)',...
            'FontSize',14,'FontWeight','bold')
    end
    if cl >= 4 && cl<=6; 
        xlabel('Age (months)','FontSize',14,'FontWeight','bold')
    end
    set(gca,'FontSize',14,'FontWeight','bold')
end
% tightfig(fig_hndl)

% ip address of baby room EEG computer
% 172.25.152.22