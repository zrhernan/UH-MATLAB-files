function classtargetvector_plot( classlabel, tasklabels )
%PLOTCLASSLABEL plots the class target vector against time samples
%   Inputs:
%       - classlabel = [n X 1] vector of time samples containing a class
%       identifier for each sample in time
%   Output:
%       - a figure plotting the target vector of classes against time
%       samples
figure;

num_classes=max(unique(classlabel));
num_samples = size(classlabel,1);

plot(classlabel,'Linewidth',2); 

% Modify plot parameters
set(gca,'YTick',[1:num_classes],...
    'YTickLabel',tasklabels,...
    'FontSize',12,'FontWeight','bold');
ylabel('Class Label','FontSize',14); ylim([0 (num_classes + 0.5)]); 
xlabel('Time Samples','FontSize',14);  xlim([0 num_samples]);

% Modify figure parameters
set(gcf, 'color', 'w');
end

