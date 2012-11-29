function [ ] = drawcircle( center, radius, linewidth )
%DRAWCIRCLE Summary of this function goes here
%   Detailed explanation goes here

if ~ exist('linewidth')
    linewidth = 1;
end

angles = [(1:5:360)*pi/180,0]';
pts = [center(1) + radius * cos(angles), center(2) + radius * sin(angles) ];
line(pts(:,1),pts(:,2), 'LineWidth', linewidth);

end


