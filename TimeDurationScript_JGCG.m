clc, clear all, close all
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
month = [9, 20, 6.1, 9.1, 6, 10, 16, 23, 18];
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
errorbar(month, TimeLengthData.meanDurationsbyInfant/60000,...
    TimeLengthData.stdDurationsbyInfant/60000,'.r',...
    'markersize', 24, 'linewidth', 2)
ylabel('Average Behavior Duration (minutes)',...
            'FontSize',14,'FontWeight','bold')
xlabel('Age (months)','FontSize',14,'FontWeight','bold')
set(gca,'FontSize',14,'FontWeight','bold')


%{
figure
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
%}

%% boxplots per infant

for i = 1:length(InfantID)
    B{i} = TimeLengthData.durationsByInfant{i};
end

for i = 1:length(B)
    len(i) = length(B{i});
end
L = max(len);

Box = nan(L,length(B));

for i = 1:length(B)
    beh = B{i}./1000;
    Box(1:len(i),i) = beh;
end
figure
boxplot(Box,month,'plotstyle','compact')
ylabel('behavior duration (s)')
xlabel('subject')



%% box plots by behavior by infant
labelss = {'6','6','9','9','10','16','18','20','23'};

for i = 1:length(InfantID)
    for j = 1:6
        BI{i,j} = TimeLengthData.durations{i,j};
    end
end


for i = 1:length(BI)
    for j = 1:6
        lenI(i,j) = length(BI{i,j});
    end
end
L = max(max(lenI));

BoxI = nan(L,length(BI),6);

for i = 1:length(BI)
    for j = 1:6
        beh = BI{i,j}./1000;
        BoxI(1:lenI(i,j),i,j) = beh;
    end
end

figure
for i = 1:6
    subplot(1,6,i)
    boxplot(BoxI(:,:,i),month,'plotstyle','compact','symbol','');
    ylabel('behavior duration (s)')
    %xlabel('subject')
    title([TimeLengthData(1).BehaviorName{i}])
    axis([0.5 9.5 0 25])
    set(gca,'Xticklabel',labelss)
end

for i = 1:6
    for j = 1:9
        numm(i,j) = length(BoxI(:,j,i))-sum(isnan(BoxI(:,j,i)));
    end
end

A = [month;numm];
out = sortrows(A',1)';

figure
ha = tight_subplot(1,6,[.01 .01],[.1 .1],[.04 .01]);
for i = 1:6
    % subplot(1,6,i)
    %    ha = tight_subplot(1,6,[.01 .03],[.1 .01],[.01 .01]);
    axes(ha(i));
    bar(1:9,out(i+1,:));
    if i ==1
       ylabel('T behavior duration (s)')
       else
       ylabel('')
       set(gca,'yticklabel','')
    end
    title([TimeLengthData(1).BehaviorName{i}])
    axis([0.5 9.5 0 75])
    set(gca,'Xticklabel','')
end


figure
for i = 1:6
   subplot(1,6,i)
    bar(1:9,out(i+1,:));
    if i ==1
    ylabel('behavior duration (s)')
    else
    ylabel('')
    set(gca,'yticklabel','')
    end
   % xlabel('subject')
    title([TimeLengthData(1).BehaviorName{i}])
   axis([0.5 9.5 0 75])
   set(gca,'Xticklabel',labelss)
end




%
% Plot them together
figure
ha = tight_subplot(2,6,[.01 .01],[.01 .01],[.04 .01]);
for i = 1:6
    % subplot(1,6,i)
    %    ha = tight_subplot(1,6,[.01 .03],[.1 .01],[.01 .01]);
    axes(ha(i));
    bar(1:9,out(i+1,:));
    if i ==1
       ylabel('T behavior duration (s)')
       else
       ylabel('')
       set(gca,'yticklabel','')
    end
    title([TimeLengthData(1).BehaviorName{i}])
    axis([0.5 9.5 0 75])
    set(gca,'Xticklabel','')
end

ha1 = tight_subplot(2,6,[.02 .01],[.01 .01],[.04 .01]);
for i = 7:12
    axes(ha1(i));
     %bar(1:9,out(i+1-6,:));
    boxplot(BoxI(:,:,i-6),month,'plotstyle','compact','symbol','');
    if i == 7
        ylabel('behavior duration (s)')
    else
        ylabel('')
        set(gca,'yticklabel','')
    end
    %xlabel('subject')
    %title([TimeLengthData(1).BehaviorName{i-6}])
    axis([0.5 9.5 0 25])
    set(gca,'Xticklabel','')
end



%
% Plot them together
figure
for i = 1:6
    subplot(1,6,i)
    bar(1:9,out(i+1,:),'k');
    if i ==1
       ylabel('Number of events')
       else
       ylabel('')
       set(gca,'yticklabel','')
    end
    title([TimeLengthData(1).BehaviorName{i}])
    axis([0.5 9.5 0 75])
    set(gca,'Xticklabel','')
end

figure
for i = 1:6
    subplot(1,6,i)
    boxplot(BoxI(:,:,i),month,'symbol','','color','k');
    if i == 1
        ylabel('Duration (s)')
        %set(gca,'yticklabel','')
    else
        ylabel('')
        set(gca,'yticklabel','')
    end
    %xlabel('subject')
    %title([TimeLengthData(1).BehaviorName{i-6}])
    %axis([0.5 9.5 0 25])
    set(gca,'Xticklabel',labelss)
end



%% Imitate figure

figure
subplot(2,1,1)
bar(1:9,out(3,:),'k');
set(gca,'Xticklabel','')
ylim([0 50])
xlim([0.5 9.5])
ylabel('Number of events')
subplot(2,1,2)
boxplot(BoxI(:,:,2),month,'symbol','','color','k');
ylabel('Duration (s)')
ylim([0 20])
xlim([0.5 9.5])
text(4,50,'Behavior: Imitate') 



%% Pie charts
BoxIp = BoxI;
[AA,ind] = sort(month);
for i = 1:size(BoxIp,1)
    for j = 1:size(BoxIp,2)
        for k = 1:size(BoxIp,3)
            if isnan(BoxIp(i,j,k))
                BoxIp(i,j,k) = 0.0000001;
            end
        end
    end
end


indice = [2 3 6 7 8 9];

labs = {'','','','','',''};

neworder = [2 5 6 1 3 4];
legs = (TimeLengthData(1).BehaviorName(neworder));
explode=[1 0 0 0 0 0];
figure
for i = 1:9
    subplot(3,3,i)
    xx = nansum(BoxIp(:,ind(i),:));
    X(:,i) = squeeze(xx./sum(xx));
    XX(:,i) = X(neworder,i);
    pie(XX(:,i),explode,labs);
    %legend(TimeLengthData(1).BehaviorName,'Location','southoutside','Orientation','horizontal')
    
      if i == 1
          legend(legs)
      end
    
    title(InfantID(ind(i)))
end










