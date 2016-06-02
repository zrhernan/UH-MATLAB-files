function [ data1,data2,C1,C2 ] = kFoldXvalidation( XvalidList,iteration,fmrx_1,fmrx_2,fmrx_3,ic_1,ic_2,ic_3)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
% Specify 6 equal partitions of trial block data
fmrx_1hf1=fmrx_1(1:round((50/100)*length(fmrx_1(:,1))),:);
fmrx_1hf2=fmrx_1(round((50/100)*length(fmrx_1(:,1))):end,:);
fmrx_2hf1=fmrx_2(1:round((50/100)*length(fmrx_2(:,1))),:);
fmrx_2hf2=fmrx_2(round((50/100)*length(fmrx_2(:,1))):end,:);
fmrx_3hf1=fmrx_3(1:round((50/100)*length(fmrx_3(:,1))),:);
fmrx_3hf2=fmrx_3(round((50/100)*length(fmrx_3(:,1))):end,:);
ic_1hf1=ic_1(1:round((50/100)*length(ic_1(:))));
ic_1hf2=ic_1(round((50/100)*length(ic_1(:))):end);
ic_2hf1=ic_2(1:round((50/100)*length(ic_2(:))));
ic_2hf2=ic_2(round((50/100)*length(ic_2(:))):end);
ic_3hf1=ic_3(1:round((50/100)*length(ic_3(:))));
ic_3hf2=ic_3(round((50/100)*length(ic_3(:))):end);
 % Training data
    eval(['train_fmrx=[fmrx_',num2str(XvalidList(iteration,1)),';fmrx_',num2str(XvalidList(iteration,2)),'];']);         %,';fmrx_',num2str(XvalidList(iteration,3)),'];']);
    eval(['train_ic=[ic_',num2str(XvalidList(iteration,1)),';ic_',num2str(XvalidList(iteration,2)),'];']);          %,';ic_',num2str(XvalidList(iteration,3)),'];'])
    [data2,C2] = ytoc(train_fmrx',train_ic);
  % Testing data
    eval(['[data1,C1] = ytoc(transpose(fmrx_',num2str(XvalidList(iteration,3)),'),ic_',num2str(XvalidList(iteration,3)),');']);

end

