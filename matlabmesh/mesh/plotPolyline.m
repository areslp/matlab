function [ ] = plotPolyline( P, indices )
%[ ] = plotPolyline( P, indices )
%   Detailed explanation goes here

borderWidth = 3;
boundaryColor = 'b';

vi = [indices(:); indices(1)];
if size(P,2) == 3
    line( P(vi,1), P(vi,2), P(vi,3), 'LineWidth', borderWidth, 'Color', boundaryColor );
else
    line( P(vi,1), P(vi,2), 'LineWidth', borderWidth, 'Color', boundaryColor );
end
%scatter3( P(vi,1), P(vi,2), P(vi,3), 'r', 'filled');

end
