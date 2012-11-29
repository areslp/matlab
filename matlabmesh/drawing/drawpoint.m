function [ ] = drawpoint( pt, pointsize )
%DRAWCIRCLE draw a point 
%   Detailed explanation goes here

if ~ exist('pointsize')
    pointsize = 0.01;
end
    

rectangle('Position',[pt(1)-pointsize/2,pt(2)-pointsize/2,pointsize,pointsize],'Curvature',[1,1]);

