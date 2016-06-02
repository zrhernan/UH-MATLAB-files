clc; clear all; close all;

%% Annotation for Infant Data
InfantDataAnnotList = {'A06-09-28-2013','B06-10-30-2013','N09-12-04-2013'...
    'GR09-07-12-2014','LW10-06-19-2014','AR16-07-14-2014','J20-08-19-2013','RB23-12-04-2014'};

%% Importation of Class Labeling Vector
cd('\\bmi-nas-01\Contreras-UH\Infantdata\Infantdata\code\Zachs_Infant_decoding_files')
for infant = 1:8;
    load(['Data\',InfantDataAnnotList{infant},'\classlabel.mat']);     % load class vector
    TargetVectors{infant,1} = classlabel;
    TargetVectors{infant,2} = ((length(find(classlabel~=0)))/100);
end

%% change file path for EEG usage later
%{
cd('C:\Users\jgcruz\Documents\MATLAB\eeglab13_4_4b')
eeglab
%}

%% Histogram Age
A = [6;7;18;16;6;13;17;9;20;10;24;18;10;16;9;21;22;9;23;6;15;8;7; ...
    17;6;6;21;9;14;10;13;10;6;10;18;16;7;15;8;24;21;20;14;7;8;22];


figure; 
X = A; 
nbins = 24; 
hist(X,[6:1:24]); 
axis([5.5 24.5 0 8]); 
fontsize = 14; 
legend('n=47'); 
xlabel('Age (months)','FontSize',fontsize); 
ylabel('Number of Subjects','FontSize',fontsize); 
title('Histogram Age of Infants recruited','FontSize',fontsize); 
set(gca,'FontSize',fontsize,'Xtick',2:2:24); 

%% Histogram Gender
M = 23;
F = 24;
figure; bar([M F]); 
axis([-.2 3.2 0 30]); 
xlabel('Gender','FontSize',fontsize); 
ylabel('Number of Subjects','FontSize',fontsize); 
legend('n=47'); 
set(gca,'FontSize',fontsize); 
title('Histogram Gender'); 
set(gca,'XTickLabel',{'M','F'})


%% Behavior segmentation histograms

% 8 subjects
B = ...
[[19	34	36	42	18	32	36	65];
[0	0	4	0	1	23	24	17];
[17	17	23	46	17	39	28	73];
[12	5	30	21	10	8	36	21];
[0	0	14	2	2	27	17	46];
[2	6	11	4	1	0	5	16];
[50	62	118	115	49	129	146	238]];  % B = behavior x subject

Etimes = [1460.81; 1320.478; 958.31; 798.75; 580.76; 946.47; 895.35; 1504.02];



Subjectlabels = {'A06 (6m)','B06 (6m)','N09 (9m)','GR09 (9m)','LW10 (10m)','AR16 (16m)','J20 (20m)','RB23 (23m)'};
SubjectAge = [06 06 09 09 10 16 20 23];
Behaviorlabels = {'Reach-grasp','Reach-offer','Observe','Explore','Imitate','Rest'};

figure;
for i = 1:8
    subplot(2,4,i)
    bar(B(1:end-1,i))
    title(Subjectlabels(i))
    set(gca,'xticklabel',Behaviorlabels,'fontsize',10)
    h=gca;
    th = rotateticklabel(h,30);
    legend(['n=',num2str(B(end,i))]);
    axis([0 7 0 1.20*max(B(1:end-1,i))])
end

%% MANOVA
y1 = B(end,:)';
y2 = Etimes;
x1 = SubjectAge;

[d,p] = manova1([y1 y2],x1);

%% Number of segmented behaviors per baby's age (before he gets tired)

fts= 20;

figure
scatter(SubjectAge,B(end,:),70,'filled','k')
axis([5.1 24.9 0 300])
xlabel('Age (months)','fontsize',fts)
ylabel('Number of behaviors segmented','fontsize',fts,'fontweight','bold')
title('Number of identified behaviors vs infant''s age','fontsize',fts,'fontweight','bold')
set(gca,'fontsize',fts,'fontweight','bold')
set(gca,'xtick',6:2:24)
hold on

x = SubjectAge;
y = B(end,:);
p = polyfit(x,y,1);
yfit = polyval(p,x);
yresid = y - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(y)-1) * var(y);
rsq = 1 - SSresid/SStotal;
plot(x,yfit,'r')
%legend(['r^2=',num2str(rsq)])
annotation('textbox', [0.7,0.8,0.1,0.1],...
           'String', ['r^2=',num2str(rsq)],'fontsize',fts);
hold off


%% Total time of active involvement
for i = 1:size(TargetVectors,1)
    AT(i) = (TargetVectors{i,2});
end


figure
scatter(SubjectAge,AT,'filled','k')
axis([5.1 24.9 0 700])
xlabel('Age (months)','fontsize',12)
ylabel('Total time of active involvement (s)','fontsize',12)
title('Total time of active involvement vs infant''s age','fontsize',12)
set(gca,'fontsize',12)
set(gca,'xtick',6:2:24)
hold on

x = SubjectAge;
y = AT;
p = polyfit(x,y,1);
yfit = polyval(p,x);
yresid = y - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(y)-1) * var(y);
rsq = 1 - SSresid/SStotal;
plot(x,yfit,'r')
%legend(['r^2=',num2str(rsq)])
annotation('textbox', [0.7,0.8,0.1,0.1],...
           'String', ['r^2=',num2str(rsq)]);
hold off


%%
figure;
scatter(SubjectAge,B(end,:)./AT,'filled','k');
axis([5.1 24.9 0 0.6])
xlabel('Age (months)','fontsize',12)
ylabel('Behavior change rate (1/s)','fontsize',12)
title('Behavior change rate vs infant''s age','fontsize',12)
set(gca,'fontsize',12)
set(gca,'xtick',6:2:24)

%% Plot number of actions and time
figure
scatter(y1,y2);

%% Histogram of Behaviors by age group
agegroup = [[6 8];
            [9 11];
            [12 15];
            [16 20];
            [21 24]];
        
figure;
A = SubjectAge;
for i = 1:size(agegroup,1)
    L = agegroup(i,1);
    U = agegroup(i,2);
    idx{i} = [find(A>=L & A<=U)];
end

for i = 1:6
    for j = 1:size(agegroup,1)
        A = B(1:end-1,idx{j})';
        S(i,j,:) = sum(A);
    end
end




