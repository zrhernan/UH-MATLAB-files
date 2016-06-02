%% Complexity measure using Elliptic Fourier Analysis (EFA) program
% Sho Nakagome 10th Sep 2015
% Feel free to give any feedbacks: snakagome@uh.edu

clearvars -except bbs cbl dar wbg

%% Move to the folder (where data are)
cd('/Users/sho/Dropbox/UH/Research/Projects/Avater_Project/ICVR/MATLAB');

%% Load color_code.m for colors
color_code

%% Static Variables
Fs = 100; % sampling frequency 100Hz
d_pre = 0; % initialize distance
d_exp = 0;
dis_pre = 0;
dis_exp = 0;
SSE_max_p = 0;
SSE_p = 0;
SSE_max_e = 0;
SSE_e = 0;
sp = 0;
se = 0;
ee = 0;

%% Loading data
sg01_01 = load('Data/SG01-T01-15-09-29-kin.mat');
sg01_02 = load('Data/SG01-T02-15-10-27-kin.mat');
sg01_03 = load('Data/SG01-T03-15-11-05-kin.mat');
sg02_01 = load('Data/SG02-T01-15-11-02-kin.mat');
sg02_02 = load('Data/SG02-T02-15-11-06-kin.mat');
sg02_03 = load('Data/SG02-T03-15-11-11-kin.mat');
sg03_01 = load('Data/SG03-T01-15-11-02-kin.mat');
sg03_02 = load('Data/SG03-T02-15-11-06-kin.mat');
sg03_03 = load('Data/SG03-T03-15-11-11-kin.mat');

% change these variables (e.g. assign subject ID to var: subject)
subject = sg03_03;
session = 3; % which session was it?

gcycle = subject.var.conductor.tGaitCycleHeel;
gp = subject.var.kinematics.rheelpos.movdot(:,1:2);

%% Dividing the data using gait cycles
if session == 1
    % find gait landmarks - finding start of pre
    for i = 1:length(gcycle)
        if gcycle(i,1) > 420
            sp = i - 1;
            break;
        end
    end

    % find gait landmarks - finding start of exp
    for i = 1:length(gcycle)
        if gcycle(i,1) > 900
            se = i - 1;
            break;
        end
    end

    % find gait landmarks - finding end of exp
    for i = 1:length(gcycle)
        if gcycle(i,1) > 1380
            ee = i - 1;
            break;
        else
            ee = length(gcycle);
        end
    end
end

if (session == 2) || (session == 3)
    % find gait landmarks - finding start of exp
    for i = 1:length(gcycle)
        if gcycle(i,1) > 420
            se = i - 1;
            break;
        end
    end

    % find gait landmarks - finding end of exp
    for i = 1:length(gcycle)
        if gcycle(i,1) > 900
            ee = i - 1;
            break;
        else
            ee = length(gcycle);
        end
    end
end

% number of gait cycles in each phase (pre- and exp-)
if session == 1
    ng_pre = (se - 1) - sp;
end
ng_exp = ee - se;

% calculating centroids in pre-exposure phase
if session == 1
    for i = 1:ng_pre
        % calculating mean for each gait
        m_pre(i,1) = mean(gp(round(gcycle(sp+i,1)*Fs):round(gcycle(sp+i,2)*Fs),1));
        m_pre(i,2) = mean(gp(round(gcycle(sp+i,1)*Fs):round(gcycle(sp+i,2)*Fs),2));
    end
end

% doing the same for early-exposure phase
for i = 1:round(ng_exp/2)
    m_exp(i,1) = mean(gp(round(gcycle(se+i,1)*Fs):round(gcycle(se+i,2)*Fs),1));
    m_exp(i,2) = mean(gp(round(gcycle(se+i,1)*Fs):round(gcycle(se+i,2)*Fs),2));
end

% doing the same for post-exposure phase
for i = 1:round(ng_exp/2)
    m_exp_pt(i,1) = mean(gp(round(gcycle(ee+i-1-round(ng_exp/2),1)*Fs):round(gcycle(ee+i-1-round(ng_exp/2),2)*Fs),1));
    m_exp_pt(i,2) = mean(gp(round(gcycle(ee+i-1-round(ng_exp/2),1)*Fs):round(gcycle(ee+i-1-round(ng_exp/2),2)*Fs),2));
end

%% Scatter histogram
fsz_s = 15;
fsz_l = 18;
m_pre = randn(33,2)+10;
m_exp = randn(33,2)+5;
m_exp_pt = randn(33,2)+8;

session=1;
figure
% enlarge figure to full screen.
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
subtightplot(3,3,[4,5,7,8]);
if session == 1
    scatter(m_pre(:,1),m_pre(:,2),'MarkerEdgeColor',dar(2,:),...
        'MarkerFaceColor',dar(2,:)); hold on;
end
scatter(m_exp(:,1),m_exp(:,2),'*','MarkerEdgeColor',cbl(2,:),...
    'MarkerFaceColor',cbl(2,:)); hold on;
scatter(m_exp_pt(:,1),m_exp_pt(:,2),'x','MarkerEdgeColor',bbs(2,:),...
    'MarkerFaceColor',bbs(2,:)); hold on;
% legend properties
if session == 1
    cen_leg = legend('pre centroid','early exp','late exp');
elseif (session == 2) || (session == 3)
    cen_leg = legend('early exp','late exp');
end
legend('boxoff');
set(cen_leg,'FontSize',fsz_l);
% label properties
xlabel('relative x position (mm)','FontSize',fsz_l);
ylabel('relative y position (mm)','FontSize',fsz_l);
set(gca,'FontSize',fsz_s);
xl = xlim;
yl = ylim;

% properties for histogram
% transparency
tp = 0.5;
nbins_h1 = 5;
nbins_h2 = nbins_h1/((xl(2) - xl(1))/(yl(2) - yl(1)));

subtightplot(3,3,[1,2]);
if session == 1
    h_pre = histogram(m_pre(:,1));
    h_pre.FaceColor = dar(2,:);
    h_pre.FaceAlpha = tp;
    hold on;
end
h_exp = histogram(m_exp(:,1));
h_exp.FaceColor = cbl(2,:);
h_exp.FaceAlpha = tp;
hold on;
h_exp_pt = histogram(m_exp_pt(:,1));
h_exp_pt.FaceColor = bbs(2,:);
h_exp_pt.FaceAlpha = tp;
if session == 1
    h_pre.Normalization = 'probability';
    h_pre.BinWidth = nbins_h1;
end
h_exp.Normalization = 'probability';
h_exp.BinWidth = nbins_h1;
h_exp_pt.Normalization = 'probability';
h_exp_pt.BinWidth = nbins_h1;
% deleting xtick
set(gca,'XTickLabel','');
% replacing yticks with percentage
yTickSet = get(gca,'ytick');
new_yticks = [cellstr(num2str(get(gca,'ytick')'*100))];
set(gca,'YTickLabel',new_yticks);
yTickLabelSet = get(gca,'YTickLabel');
yTickLabelSet{1} = '';
set(gca,'yTick',yTickSet,'YTickLabel',yTickLabelSet);
ylabel('Percentage (%)','FontSize',fsz_l);
set(gca,'FontSize',fsz_s);
box off;

subtightplot(3,3,[6,9]);
if session == 1
    h2_pre = histogram(m_pre(:,2),'Orientation','horizontal');
    h2_pre.FaceColor = dar(2,:);
    h2_pre.FaceAlpha = tp;
    hold on;
end
h2_exp = histogram(m_exp(:,2),'Orientation','horizontal');
h2_exp.FaceColor = cbl(2,:);
h2_exp.FaceAlpha = tp;
hold on;
h2_exp_pt = histogram(m_exp_pt(:,2),'Orientation','horizontal');
h2_exp_pt.FaceColor = bbs(2,:);
h2_exp_pt.FaceAlpha = tp;
if session == 1
    h2_pre.Normalization = 'probability';
    h2_pre.BinWidth = nbins_h2;
end
h2_exp.Normalization = 'probability';
h2_exp.BinWidth = nbins_h2;
h2_exp_pt.Normalization = 'probability';
h2_exp_pt.BinWidth = nbins_h2;
set(gca,'YTickLabel','');
xTickSet = get(gca,'xtick');
new_xticks = [cellstr(num2str(get(gca,'xtick')'*100))]; 
set(gca,'XTickLabel',new_xticks);
xTickLabelSet = get(gca,'XTickLabel');
xTickLabelSet{1} = '';
set(gca,'xTick',xTickSet,'XTickLabel',xTickLabelSet);
xlabel('Percentage (%)','FontSize',fsz_l);
set(gca,'FontSize',fsz_s);
box off;