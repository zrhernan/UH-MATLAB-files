function enhanceFigure( figure_handle )
%ENHANCEFIGURE will maximize the figure to full screen, tighten the plot
%(or subplots), and re-maximize the screen again for saving the figure
%   The purpose is to provide an automated way of saving large figures
%   This uses the TIGHTFIG function created by Richard Crozier, which you
%   can find using this link:
%   {http://www.mathworks.com/matlabcentral/fileexchange/34055-tightfig}
warning off MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame
jFrame = get(handle(figure_handle),'JavaFrame'); % using Javaframe to maximize figure
pause(10)  % wait 10 seconds so JavaScript function has ample time to execute
jFrame.setMaximized(true);   % to maximize the figure
pause(10)  % wait 10 seconds so JavaScript function has ample time to execute
tightfig;  % to tighten the margins around the figure
pause(10)  % wait 10 seconds so tightfig function has ample time to execute
jFrame.setMaximized(true);   % to re-maximize the figure
pause(10)  % wait 10 seconds so JavaScript function has ample time to execute
end