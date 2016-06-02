function timecompute( toc )
%TIMECOMPUTE displays the time of a computation process using 'toc' from
%MATLAB. Displayed time given in seconds, minutes, and hours.

if toc >= 3600; % if numbers of seconds longer than an hour
    tocHours = toc/3600;
    tocMins = (tocHours - floor(tocHours))*60;
    tocSecs = (tocMins - floor(tocMins))*60;
disp(['Computation Time:: ',...
    num2str(floor(tocHours)), ' hours | ',...
    num2str(floor(tocMins)), ' minutes | '...
    num2str(floor(tocSecs)) ' seconds'])
elseif toc >= 60 && toc < 3600; % if numbers of seconds longer than a minute but less than an hour
    tocMins = toc/60;
    tocSecs = (tocMins - floor(tocMins))*60;
disp(['Computation Time:: ',...
    num2str(floor(tocMins)), ' minutes | '...
    num2str(floor(tocSecs)) ' seconds'])
elseif toc < 60; % if numbers of seconds longer than a minute but less than an hour
    tocSecs = toc;
    disp(['Computation Time:: ',...
    num2str(floor(tocSecs)) ' seconds'])    
end

