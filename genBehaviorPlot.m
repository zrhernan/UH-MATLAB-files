% Histogram generation for Ages, Gender, behaviors in a comprehensive,
% cover-it-all figure

clc; clear all; close all;

%% load data
% gender.mat, ages.mat, behaviors.mat, infantlabels.mat
ilabels = load('infantlabels');

datacol0 = load('wasdatacollected');
datacol = cell2mat(datacol0.wasdatacollected(:,2));

% Age histogram data
AgeS = load('ages');
Age = AgeS.ages;
Age1 = Age(find(datacol==1));   %number of subjects with successful data collection
Age2 = Age(find(datacol==0));   %number of subjects with no data collection


Age1M = [2 2 1 3 4]; Age2M= [1 0 1 0 1];
Age1F =[2 3 1 3 2]; Age2F=[0 0 0 0 1];

% Gender histogram data
GenderS = load('gender'); 
Genderd = GenderS.gender;
Gender = cell2mat(Genderd(:,2)); 

% Behavior histogram data
Behavior = load('behaviors'); % B = behavior x subject
Behavior = Behavior.behaviors(1:end-1,:);
Behaviorlabels = {'Reach-grasp','Reach-offer','Observe','Explore','Imitate','Rest'};

%% plot
agegroup1 = {'6-8', '9-11','12-15', '16-19', '20-24'};
nbinsag1 = [7 10 13.5 17.5 22];
agegroup2 = {'6-8', '9-11','12-14', '15-17', '18-20','21-24'};
nbinsag2 = [7 10 13 16 19 22.5];


% Age histogram

% age groups
X = Age1; 
Y = Age2;
nbins = nbinsag1; 
[p1,q1]=hist(X,nbins); 
[p2,q2]=hist(Y,nbins); 
plotthis = [p1;p2];
xthis = [q1];

figure(1); 
bar(plotthis','stacked')
axis([0.5 5.5 0 10]); 
fontsize = 14; 
legend(['n=',num2str(length(Age1))],['n=',num2str(length(Age2))]); 
%xlabel('Age (months)','FontSize',fontsize); 
ylabel('Number of Subjects','FontSize',fontsize); 
title('Histogram Age of Infants recruited','FontSize',fontsize); 
set(gca,'FontSize',fontsize,'Xtick',(1:6),'XTickLabel',agegroup1); 



% monthly
X = Age1; 
Y = Age2;
nbins = [1:24]; 
[p11,q11]=hist(X,nbins); 
[p12,q12]=hist(Y,nbins); 
plotthis2 = [p11;p12];
xthis = [q11];

figure(2); 
X = Age1; 
Y = Age2;
bar(plotthis2','stacked')
axis([5.5 24.5 0 4]); 
fontsize = 14; 
legend(['n=',num2str(length(Age1))],['n=',num2str(length(Age2))]); 
xlabel('Age (months)','FontSize',fontsize); 
ylabel('Number of Subjects','FontSize',fontsize); 
title('Histogram Age of Infants recruited','FontSize',fontsize); 
set(gca,'FontSize',fontsize,'Xtick',[1:24]); 



% Gender histogram
Mus = 3; Fus = 1; %number of subjects with unsuccessful data collection
M1 = Gender(1)-Mus; F1 = Gender(2)-Fus; %number of subjects with successful data collection
plotthisG = [[M1 F1];[Mus Fus]]';

figure(3); 
bar(plotthisG,'stacked'); 
axis([0 3 0 20]); 
xlabel('Gender','FontSize',fontsize); 
ylabel('Number of Subjects','FontSize',fontsize); 
legend(['n=',num2str(M1+F1)],['n=',num2str(Mus+Fus)]); 
set(gca,'FontSize',fontsize); 
title('Histogram Gender'); 
set(gca,'XTickLabel',{'M','F'})



% Behavior histogram
bsum = sum(Behavior,2);
avgbsum = bsum/[1];

figure(4); 
bar(bsum); 
axis([0 7 0 max(bsum)*1.1]); 
ylabel('Number of Events','FontSize',fontsize); 
title('Histogram Behaviors'); 
legend(['n=',num2str(sum(bsum))]); 
fontsize = 12;
set(gca,'FontSize',fontsize); 
set(gca,'xTicklabel',Behaviorlabels,'fontsize',10)
% h=gca;
% th = rotateticklabel(h,30);
view([-90 90]);



% Data density plot
%Calculating Behaviors by Age group
ag1 = [[6 8];[9 11];[12 15];[16 19];[20 24]];
for i = 1:length(ag1)
    idx = find(Age>=ag1(i,1) & Age<=ag1(i,2));
    bsumD(:,i) = sum(Behavior(:,idx),2);
end

 figure(5)
 imagesc(flipud(bsumD)); %view([0 90])
 hcb = colorbar;
 set(hcb,'yTick',[20:20:160],'yticklabel',{'20','40','60','80','100','120','140','Events'},'fontsize',fontsize)
 set(gca,'YTickLabel','');
 set(gca,'XTickLabel','');
 set(gca,'fontsize',14);
 box off;
 colormap winter
 

 
 %% All together
 figure
 fontsize = 12;
 
 
 % Big consolidated figure
subplot(3,3,2:3);
bsum2 = sum(bsumD,1);
H = bar(bsum2,'BarWidth', 1);
axis([0.5 5.5 0 750]); 
%legend(['n=',num2str(sum(bsum2))]); 
ylabel('Number of Subjects','FontSize',fontsize,'fontweight','bold'); 
title('Event count by age group','FontSize',fontsize,'fontweight','bold'); 
set(gca,'FontSize',fontsize,'Xtick',(1:6),'XTickLabel',''); 
x=[1:5]';
y=abs([3 5 1 4 5]);
for i1=1:numel(y)
    text(x(i1),y(i1),['N=',num2str(y(i1))],'fontsize',10,'color','w',....
               'HorizontalAlignment','center',...
               'VerticalAlignment','bottom')
end
fontsize = 12;


subplot(3,3,[4,7]);
H = bar(flipud(bsum),'BarWidth', 1); 
axis([0.5 6.5 0 max(bsum)*1.1]); 
ylabel('Number of Events','FontSize',fontsize,'fontweight','bold'); 
title({'Event count ', 'for each behavior'},'fontsize',12,'fontweight','bold'); 
%legend(['n=',num2str(sum(bsum))]); 
set(gca,'FontSize',fontsize); 
set(gca,'xTicklabel',Behaviorlabels(end:-1:1))
% h=gca;
% th = rotateticklabel(h,30);
view([-90 90]);


subplot(3,3,[5:6,8:9])
 imagesc((bsumD)); %view([0 -90])
 %colorbar('southoutside')
 set(gca,'YTickLabel','');
 set(gca,'XTickLabel','');
 set(gca,'fontsize',14);
 box off;
 colormap(jet)
xlabel('Age (months)','FontSize',fontsize,'fontweight','bold'); 
set(gca,'FontSize',fontsize,'Xtick',(1:6),'XTickLabel',agegroup1); 


% Demographics
figure
 subplot(1,3,3)
H = bar(plotthisG,'stacked'); 
set(H(1),'facecolor','b')
set(H(2),'facecolor',[0.7 0.7 0.7])
axis([0 3 0 20]); 
xlabel('Gender','FontSize',fontsize); 
ylabel('Number of Subjects','FontSize',fontsize); 
%legend(['successful test = ',num2str(M1+F1)],['unsuccessful test = ',num2str(Mus+Fus)]); 
set(gca,'FontSize',fontsize); 
title('Histogram Gender'); 
set(gca,'XTickLabel',{'M','F'})

 subplot(1,3,1:2)
H = bar(plotthis','stacked','BarWidth', 1);
set(H(1),'facecolor','b')
set(H(2),'facecolor',[0.7 0.7 0.7])
axis([0.5 5.5 0 10]); 
fontsize = 10; 
legend(['successful test = ',num2str(length(Age1))],['unsuccessful test = ',num2str(length(Age2))]); 
xlabel('Age (months)','FontSize',fontsize); 
ylabel('Number of Subjects','FontSize',fontsize); 
title('Histogram Age of Infants recruited','FontSize',fontsize); 
set(gca,'FontSize',fontsize,'Xtick',(1:6),'XTickLabel',''); 


% Complete demographics in single plot
figure
plothisComplete = [Age1M;Age1F;Age2M;Age2F]';
H = bar(plothisComplete,'stacked','BarWidth', 0.9); 
axis([0.5 5.5 0 10]); 
set(H(1),'facecolor','b')
set(H(2),'facecolor',[0.4 0.4 1])
set(H(3),'facecolor',[0.6 0.6 0.6])
set(H(4),'facecolor',[0.9 0.9 0.9])
legend('Male           N=27','Female','Male withrawn','Female withdrawn');
xlabel('Age (months)','FontSize',14); 
ylabel('Number of Subjects','FontSize',14); 
title('Subject recruitment','FontSize',14); 
set(gca,'FontSize',fontsize,'Xtick',(1:6),'XTickLabel',agegroup2); 

% subplot(4,4,[14,15])
% colorbar('southside')
% colormap winter


y=abs([3 5 1 4 5]);

bsumDavg = bsumD./repmat(y,6,1);
bsumavg1 = sum(bsumDavg,1);
bsumavg2 = sum(bsumDavg,2);

% Big consolidated figure v2 (average)
fontsize = 12;
figure
subplot(3,3,2:3);
bsum2 = sum(bsumavg1,1);
H = bar(bsum2,'BarWidth', 1);
axis([0.5 5.5 0 max(bsumavg1)*1.2]); 
%legend(['n=',num2str(sum(bsum2))]); 
ylabel('Number of Subjects','FontSize',fontsize,'fontweight','bold'); 
title('Average event count by age group','FontSize',fontsize,'fontweight','bold'); 
set(gca,'FontSize',fontsize,'Xtick',(1:6),'XTickLabel',''); 
x=[1:5]';
for i1=1:numel(y)
    text(x(i1),y(i1),['N=',num2str(y(i1))],'fontsize',10,'color','w',....
               'HorizontalAlignment','center',...
               'VerticalAlignment','bottom')
end
fontsize = 12;


subplot(3,3,[4,7]);
H = bar(flipud(bsumavg2),'BarWidth', 1); 
axis([0.5 6.5 0 max(bsumavg2)*1.2]); 
ylabel('Number of Events','FontSize',fontsize,'fontweight','bold'); 
title({'Average event count',' for each behavior'},'fontsize',12,'fontweight','bold'); 
%legend(['n=',num2str(sum(bsum))]); 
set(gca,'FontSize',fontsize); 
set(gca,'xTicklabel',Behaviorlabels(end:-1:1))
% h=gca;
% th = rotateticklabel(h,30);
view([-90 90]);


subplot(3,3,[5:6,8:9])
 imagesc((bsumDavg)); %view([0 -90])
 %colorbar('southoutside')
 set(gca,'YTickLabel','');
 set(gca,'XTickLabel','');
 set(gca,'fontsize',14);
 box off;
 colormap(winter)
xlabel('Age (months)','FontSize',fontsize,'fontweight','bold'); 
set(gca,'FontSize',fontsize,'Xtick',(1:6),'XTickLabel',agegroup1); 


figure
 imagesc((bsumDavg)); %view([0 90])
 hcb = colorbar;
 caxis([0,60])
 set(hcb,'yTick',[0:10:60],'yticklabel',{'','10','20','30','40','50','Events'},'fontsize',fontsize)
 set(gca,'YTickLabel','');
 set(gca,'XTickLabel','');
 set(gca,'fontsize',14);
 box off;
 colormap winter
 
 
 
 
rate = []; 
% Big consolidated figure v3 (average rate)
fontsize = 12;
figure
subplot(3,3,2:3);
bsum2 = sum(bsumavg1,1);
H = bar(bsum2,'BarWidth', 1);
axis([0.5 5.5 0 max(bsumavg1)*1.2]); 
%legend(['n=',num2str(sum(bsum2))]); 
ylabel('Number of Subjects','FontSize',fontsize,'fontweight','bold'); 
title('Average event count by age group','FontSize',fontsize,'fontweight','bold'); 
set(gca,'FontSize',fontsize,'Xtick',(1:6),'XTickLabel',''); 
x=[1:5]';
for i1=1:numel(y)
    text(x(i1),y(i1),['N=',num2str(y(i1))],'fontsize',10,'color','w',....
               'HorizontalAlignment','center',...
               'VerticalAlignment','bottom')
end
fontsize = 12;


subplot(3,3,[4,7]);
H = bar(flipud(bsumavg2),'BarWidth', 1); 
axis([0.5 6.5 0 max(bsumavg2)*1.2]); 
ylabel('Number of Events','FontSize',fontsize,'fontweight','bold'); 
title({'Average event count',' for each behavior'},'fontsize',12,'fontweight','bold'); 
%legend(['n=',num2str(sum(bsum))]); 
set(gca,'FontSize',fontsize); 
set(gca,'xTicklabel',Behaviorlabels(end:-1:1))
% h=gca;
% th = rotateticklabel(h,30);
view([-90 90]);


subplot(3,3,[5:6,8:9])
 imagesc((bsumDavg)); %view([0 -90])
 %colorbar('southoutside')
 set(gca,'YTickLabel','');
 set(gca,'XTickLabel','');
 set(gca,'fontsize',14);
 box off;
 colormap(winter)
xlabel('Age (months)','FontSize',fontsize,'fontweight','bold'); 
set(gca,'FontSize',fontsize,'Xtick',(1:6),'XTickLabel',agegroup1); 


figure
 imagesc((bsumDavg)); %view([0 90])
 hcb = colorbar;
 caxis([0,60])
 set(hcb,'yTick',[0:10:60],'yticklabel',{'','10','20','30','40','50','Events'},'fontsize',fontsize)
 set(gca,'YTickLabel','');
 set(gca,'XTickLabel','');
 set(gca,'fontsize',14);
 box off;
 colormap winter
 
%{
%% Sample code

% This is just counting up totals for gender, ages, and preferences of piece based on the %variables specific to our study
survey = cell2mat(QuestionnaireDataFINAL(2:end,12:14));
subject = QuestionnaireDataFINAL(2:end,1);
SEX = LogSheetFINAL(2:end,4);
male = strcmp(SEX,'M'); SEX_temp(male) = 1;female = strcmp(SEX,'F');SEX_temp(female) = 2; SEX_temp = SEX_temp';
clear SEX male female
for i = 1:length(subject);
    visual(i) = survey(i,2);
    emotion(i) = survey(i,3);
    idx = find(strcmp(LogSheetFINAL(:,1),subject{i}));
    age(i) = LogSheetFINAL{idx,9};
    sex(i) = SEX_temp(idx-1);
end
for i = 1:8; visplease_bar(i) = length(find(survey(:,2) == i)); emotion_bar(i) = length(find(survey(:,3) == i));end
male = length(find(SEX_temp == 1));female = length(find(SEX_temp == 2));


% This is where the above data is used in the plot.
figure;
subplot(3,4,2:3);
hist(age,50);set(gca,'fontsize',14);
xlim([min(age) max(age)]);
title('Age histogram','FontSize',14)
ylabel({'Num. of','participants'},'FontSize',14);
set(gca,'XTickLabel','');

%%%%%%%%%%%%%%%%%%%
subplot(3,4,[5,9]);

for i = 1:8;ind = find(emotion == i);temp = sex(ind);plotthis(i,1) = length(find(temp == 1));plotthis(i,2) = length(find(temp == 2));end
bar(plotthis,0.5,'stacked');view([-90 90]);%legend('Male','Female');
title({'Which piece did you find';'the most emotionally stimulating?'},'FontSize',14)
ylabel('Num. of participants','FontSize',14)
set(gca,'fontsize',14)
set(gca,'XTick',[1 2 3 4 5 6 7 8]);
xlim([0 8.5])

%%%%%%%%%%%%
 subplot(3,4,[6:7,10:11])

 DataDensityPlot(age,emotion,10);
 xlabel('Age (years)','FontSize',14);
 %colorbar('southoutside')
 set(gca,'YTickLabel','');
 set(gca,'fontsize',14);
 box off;ylim([-5 256])
%[],'fontsize',14,'YAxisLocation','right');


%%%%%%%%%
 subplot(3,4,[8,12])

for i = 1:8;ind = find(visual == i);temp = sex(ind);plotthis(i,1) = length(find(temp == 1));plotthis(i,2) = length(find(temp == 2));end
bar(flipud(plotthis),0.5,'stacked');view([90 90]);%legend('Male','Female','Location','best');
ylabel('Num. of participants','FontSize',14)
title({'Which piece did you find';'the most aesthetically pleasing?'},'FontSize',14)
set(gca,'fontsize',14)
set(gca,'XTick',[1 2 3 4 5 6 7 8]);set(gca,'XTickLabel',[8:-1:1]);
xlim([0.5 9])
%}