%% Extraction of Impedance data from the Brain Vision .txt files
% clc, clear all, close all

% for calling helper functions
addpath(genpath('\\bmi-nas-01\Contreras-UH\Infantdata\Infantdata\code\Zachs_Infant_decoding_files')) 
tic;

%% Annotation for Infant Data
files = dir('\\bmi-nas-01\Contreras-UH\Infantdata\Infantdata\Data\');   % assume starting from current directory
filenames = {files.name};
subdirs = filenames([files.isdir]);
expr = '-\d{4}$'; % contains the year at the end of the folder name
cnt = 1;
for s = 1:length(subdirs)
    if ~isempty(regexp(subdirs{s},expr, 'once'))
        InfantDataAnnotList{cnt} = subdirs{s};
        hypenIDX = regexp(subdirs{s},'-');
        InfantID{cnt} = subdirs{s}(1:hypenIDX(1)-1);
        disp(['Infant folder recognized: ',InfantDataAnnotList{cnt}])
        cnt = cnt + 1;
    end
end

%% Selection of infants
% infant = menu(['Which infant would you',10,' like to analyze?'],InfantID{:});
% close
load('BrainVision_1020_64ChannelOrder.mat')
C = length(channelOrder);       N = length(InfantDataAnnotList);
Infant_Zstart = zeros(C,N);
Infant_Zend = zeros(size(Infant_Zstart));

for infant = 1:length(InfantDataAnnotList);
    disp(['Acquiring impedance data from ', InfantID{infant}])
%% Path directory initializations
    serverPath1 = ['\\172.27.216.40\Contreras-UH\Infantdata\Infantdata\Data\',...
    InfantDataAnnotList{infant},'\EEG'];

%% Define impedance file paths (both 
    cd(serverPath1)
    impedancefileStart = [serverPath1,'\',InfantDataAnnotList{infant},'-start.txt'];
    impedancefileEnd = [serverPath1,'\',InfantDataAnnotList{infant},'-end.txt'];
    impedancefilePaths = {impedancefileStart; impedancefileEnd};

%% Import Impedance Values
    [~,Zchnlist] = createhighZchnlist(impedancefilePaths); % list of impedance values in a cell format {startZ,endZ}

%% Skip subject if no impedance values exist (at start and end)
    if all(isnan(Zchnlist{1,1}) & isnan(Zchnlist{1,2}))
        InfantID{infant} = [];
        Infant_Zstart(:,infant) = Zchnlist{1};
        Infant_Zend(:,infant) = Zchnlist{2};        
        continue
    end
%% Add to the Start/End Impedance-to-Infant Array
    Infant_Zstart(:,infant) = Zchnlist{1};
    Infant_Zend(:,infant) = Zchnlist{2};

end

% take the difference between impedance values
Infant_Zdiff = Infant_Zend - Infant_Zstart;

% remove subject with no recorded impedances
SubjectLabels = InfantID;
SubjectLabels(cellfun(@(InfantID) isempty(InfantID),InfantID))=[];
Infant_Zstart(:,any(isnan(Infant_Zstart)))=[];
Infant_Zend(:,any(isnan(Infant_Zend)))=[];
Infant_Zdiff(:,any(isnan(Infant_Zdiff)))=[];
N1 = length(SubjectLabels);
age = zeros(1,N1);
for i = 1:N1
    expr1 = regexp(SubjectLabels{i},'\d{2}$');
    age(i) = str2double(SubjectLabels{i}(expr1:end));
end
[~,age_ind] = sort(age,'ascend');
SL_sorted = SubjectLabels(age_ind);
chan_toposort = [1 2 33 34 35 36 3 37 4 38 5 39 6 40 7 41 42 8 43 9 10 44,...
    11 45 46 12 47 13 48 14 49 15 50 16 17 51 18 52 19 53 20 54 21 55 22,...
    23 56 24 57 25 58 26 59 27 60 61 62 63 64 28 29 30 31 32];
chanOrd_sorted = channelOrder(chan_toposort);
Zstart_sorted = Infant_Zstart(chan_toposort,age_ind);

%% Plot the Impedance Data (Experiment Start)
figure;

subtightplot(3,3,[4,5,7,8]);

set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
imagesc(Zstart_sorted);

intrvl = floor(0.1*N1);
set(gca,'fontsize',8,'fontweight','bold','ytick',(1:C),'yticklabel',channelOrder)
ylim([0.5 C+0.5]); xlim([0.5 N1+0.5])
ylabel('EEG Channels','fontsize',20,'fontweight','bold','rotation',90)
xticklabel_rotate((1:N1),90,SL_sorted,'fontsize',8,'fontweight','bold')
xlabel('Subject ID','fontsize',20,'fontweight','bold')
load('impedance_colormap_values.mat')
caxis([0 100])
cmap = colormap(impedance_colormappings);
cb_hndl = colorbar('FontWeight','bold','FontSize',14,'Location','southoutside');
ylabel(cb_hndl,'Impedance Values (k\Omega)','fontsize',20,'fontweight','bold')
% colorbar('off')
%set(gca,'position',[0.08 0.11 0.85 0.88]);

%% Plot average in x and y axis

Zstart_sorted_cols = mean(Zstart_sorted,1);
Zstart_sorted_rows = mean(Zstart_sorted,2);


subtightplot(3,3,[1,2])
bar(Zstart_sorted_cols,1);
xlim([0.5 length(Zstart_sorted_cols)+0.5])

subtightplot(3,3,[6,9])
bar(Zstart_sorted_rows,1);
view(90, 90)
xlim([0.5 length(Zstart_sorted_rows)+0.5])
line([32.5 32.5],get(hax,'YLim'),'Color',[1 0 0])










