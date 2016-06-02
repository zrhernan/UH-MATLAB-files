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
