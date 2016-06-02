
clc, clear all, close all
tic  % start computation time

%% Annotation for Infant Data
InfantDataAnnotList = {'N09-12-04-2013','J20-08-19-2013','B06-10-30-2013',...
    'GR09-07-12-2014','A06-09-28-2013','LW10-06-19-2014','AR16-07-14-2014',...
    'RB23-12-04-2014'};
infant = menu(['Which infant would you',10,' like to analyze?'],...
    'N09','J20','B06','GR09','A06','LW10','AR16','RB23');
close

%% Initialization of Data Structures
H_INFINITYFILTER = struct('EEG',{},'KINE',{},'CLASS',{},'CLASSIFIER',{});
% PROCESS(1).KINE = struct('InfantH',{},'InfantArmL',{},'InfantArmR',{},'ActorArmL',{},'ActorArmR',{});

%% Initializing and Assigning directory paths....
InfantData_dir = 'C:\Users\zrhernan\Infant_decoding_files\';
cd(InfantData_dir)

%% Set directory paths for 'HinfFilter_Lab_Share' function
addpath('\\172.27.216.40\Contreras-UH\NeuroREX\Hinfinity_Filter_TV_Weight_Estimation')

%% Importation of EEG Data
H_INFINITYFILTER(1).EEG = load(['Data\',InfantDataAnnotList{infant},'\EEGfiles.mat']);  % load structure of EEG attributes
disp('____EEG data imported')

% Gamma = [1.15,1.1,1.05];Q = [1e-10,1e-9,1e-8,1e-7,1e-6];
% start = 1; stop = 40000;
% 
% for i = 1:length(Gamma)
%     for j = 1:length(Q)
% tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization of Input Variables for H-infinity filter
EEGDataIN = H_INFINITYFILTER.EEG.uVdata_syncd;
RefChns(1) = find(strcmp(H_INFINITYFILTER.EEG.channelOrder,'Fp1')); % to find the correct index for 'Fp1'
RefChns(2) = find(strcmp(H_INFINITYFILTER.EEG.channelOrder,'Fp2')); % to find the correct index for 'Fp1'
RefChns(3) = find(strcmp(H_INFINITYFILTER.EEG.channelOrder,'FT9')); % to find the correct index for 'Fp1'
RefChns(4) = find(strcmp(H_INFINITYFILTER.EEG.channelOrder,'FT10')); % to find the correct index for 'FT10'
EEGDataIN(:,RefChns) = [];
[num_tsamp, num_rmvdEEGchns] = size(EEGDataIN);
Ref(:,1) = EEGDataIN(:,RefChns(1)) - EEGDataIN(:,RefChns(2));
Ref(:,2) = EEGDataIN(:,RefChns(3)) - EEGDataIN(:,RefChns(4));
Ref(:,3) = ones(num_tsamp,1);

P0 = 0.5*eye(size(Ref,2));
Pt = cell(1,num_rmvdEEGchns);
for m = 1:num_rmvdEEGchns; Pt{m} = P0; end
wh = zeros(size(Ref,2), num_rmvdEEGchns);
gamma = 1.15;  %Gamma(i);%
q = 1e-9;  %Q(j);%
[shsh] = HinfFilter_Lab_Share(EEGDataIN, Ref, gamma, Pt, wh, num_rmvdEEGchns, q);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% toc
% Yf = Yf';shsh = shsh';
% figure(i); subplot(5,1,j);
plot(EEGDataIN(:,3),'k','LineWidth',2);
hold on;
plot(shsh(:,3),'r','LineWidth',2);
ylim([-100 200]);
title(['gamma = ',num2str(gamma),', q = ',num2str(q)],'FontSize',16,'FontWeight','bold')
%     end
% end
set(gca,'FontSize',12,'FontWeight','bold')
%% Print out computation time
disp(['Computation Time:: ' num2str(toc) ' seconds or '...
    10 '.......' num2str(toc/60) ' minutes or '...
    10 '.......' num2str(toc/3600) ' hours'])