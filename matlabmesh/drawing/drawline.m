function [ ] = drawline( p0, p1, linewidth )
%DRAWCIRCLE Summary of this function goes here

if ~ exist('linewidth')
    linewidth = 1;
end

points=[p0;p1];

line( points(:,1), points(:,2), 'LineWidth', linewidth );

end