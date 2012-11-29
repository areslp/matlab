function [ faceColor ] = faceColors( mesh, faceValues, vmin, vmax, colorLow, colorHigh )
%FACECOLORS Summary of this function goes here
%   Detailed explanation goes here

alpha = (faceValues-vmin) / (vmax-vmin);
omalpha = ones(size(faceValues)) - alpha;
faceColor = [omalpha*colorLow(1) + alpha*colorHigh(1), omalpha*colorLow(2) + alpha*colorHigh(2), omalpha*colorLow(3) + alpha*colorHigh(3) ];

end
