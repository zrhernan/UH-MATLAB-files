function [ NormCM ] = ConfusionMatrix_display( confusionmatrix, tasklabels )
%ConfusMatwithPrecVal takes the confusion matrix sample number and
% calculates and displays the training set classification percentage per 
% square 
%   STEP 1: Calculate a 'normalized confusion matrix that is based on a
%   percentage of training set sizes per class
%   STEP 2: Generate confusion matrix plot
%   STEP 3: Overlay those normalized values onto the confusion matrix plot

%% STEP ONE: Calculate Normalized Confusion Matrix
% Extract only needed info from the confusion matrix
num_classes = size(confusionmatrix,1)-1; % get real class size
OA = num2str((confusionmatrix(num_classes+1,num_classes+1)*100),'%0.1f'); % get the overall accuracy value
CMsamples = confusionmatrix(1:num_classes,1:num_classes); % extract the actual confusion matrix values
% Normalize by training set size per class
trainset_size = sum(sum(CMsamples),2)/num_classes;
NormCM = (CMsamples./trainset_size)*100;

%% STEP TWO: Generate the confusion matrix plot
imagesc(NormCM);figure(gcf);
axis image;
set(gca,'XTickLabel',{},'YTickLabel',{});
set(gcf,'color','w');
colormap(flipud(gray(64)));
%% STEP THREE: Add the overlay of normalized confusion values
textStrings = num2str(NormCM(:),'%0.1f');  %# Create strings from the matrix values
textStrings = strtrim(cellstr(textStrings));  %# Remove any space padding
for txtstr = 1:size(textStrings); 
    textStrings{txtstr} = [textStrings{txtstr} '%']; % add the "%" symbol at end of each percentile value
end 
[x,y] = meshgrid(1:num_classes);   %# Create x and y coordinates for the strings
hStrings = text(x(:),y(:),textStrings(:),...      %# Plot the strings
'HorizontalAlignment','center');
midValue = mean(get(gca,'CLim'));  %# Get the middle value of the color range
%  Choose white or black for the text color of the strings so they can be 
%   easily seen over the background color
textColors = repmat(NormCM(:) > midValue,1,3);
% Change the text such it's easier to view
set(hStrings,{'Color'},num2cell(textColors,2),...
    'FontSize',24,...
    'FontWeight','bold');  

% Add text for % overall accuracy at end. Can move around to look nice on
% confusion matrix plot
annotation('textbox',...
    [0.7223 0.0596 0.0496 0.0503],...
    'String',[OA '%'],...
    'FontWeight','bold',...
    'FontSize',24,...
    'FitBoxToText','off',...
    'LineStyle','none');

%% STEP FOUR: Add Labels for the confusion matrix
%{
for k = 1:num_classes;
    text(x(k)-0.75,y(k),tasklabels(k),...      % Plot the horizontal labels
'HorizontalAlignment','right','FontSize',24,'FontWeight','bold');
end
num_classes=5;
for k = 1:num_classes;
%     if mod(k,2) == 0;
%         text(x(k)+(k-1),y(1)-0.75,tasklabels(k),...      % Plot the vertical labels
%           'HorizontalAlignment','center','FontSize',24,'FontWeight','bold','Rotation',0);  
%     else
        text(x(k)+(k-1),y(1)-0.75,tasklabels(k),...      % Plot the vertical labels
          'HorizontalAlignment','center','FontSize',24,'FontWeight','bold','Rotation',23); 
%     end
end
%}
end         %EOF

