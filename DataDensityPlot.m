function [ f ] = DataDensityPlot( x, y, levels )
%DATADENSITYPLOT Plot the data density 
%   Makes a contour map of data density
%   x, y - data x and y coordinates
%   levels - number of contours to show
%
% By Malcolm Mclean
%
    map = dataDensity(x, y, 280, 280);%256 x 256 default
    map = map - min(min(map));
    map = floor(map ./ max(max(map)) * (levels-1));
    %f = figure();
    
    image(flipud(map));
    colormap(jet(levels));
    set(gca, 'XTick', [1 40 105 170 235 280]);
    set(gca, 'XTickLabel', [min(x) 15 30 45 60 max(x)]);
    set(gca, 'YTick', [1  44 86 129 172 215 257 300]);
    set(gca, 'YTickLabel', fliplr([1 2 3 4 5 6 7 max(y)]));
    %uiwait;
end

