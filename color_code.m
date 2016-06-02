%% Color code for making figures
% 
% MATLAB R2015a
% 
% DESCRIPTION
% From: https://designschool.canva.com/blog/100-color-combinations/
% Use http://www.colorhexa.com/ to see more details for colors
%
% Copyright (c) 2015 Sho Nakagome
% snakagome@uh.edu
% All Rights Reserved
% Revision: 1.0 Date: 2015/11/04

%% %% Monochromatic color schemes
%% 03. Dark & Earthy
dar(1,:) = [.275, .129, .102];
dar(2,:) = [.412, .239, .239];
dar(3,:) = [.729, .333, .212];
dar(4,:) = [.643, .22, .125];

%% 05. Cool Blues
cbl(1,:) = [0, .231, .275];
cbl(2,:) = [.027, .341, .357];
cbl(3,:) = [.4, .647, .678];
cbl(4,:) = [.769, .875, .902];

%% 07. Watery Blue-Greens
wbg(1,:) = [.008, .11, .118];
wbg(2,:) = [0, .267, .271];
wbg(3,:) = [.173, .471, .451];
wbg(4,:) = [.435, .725, .561];

%% 30. Berry Blues
bbs(1,:) = [.118, .122, .149];
bbs(2,:) = [.157, .212, .333];
bbs(3,:) = [.302, .392, .553];
bbs(4,:) = [.816, .882, .976];

%% Plotting to see colors
% Change this variable
color = bbs;

for i = 1:length(color)
    rectangle('Position',[1+(i-1)*10,1,10,5],'FaceColor',color(i,:));
    hold on;
end