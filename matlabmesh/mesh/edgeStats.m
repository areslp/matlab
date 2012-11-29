function [ avg_len, min_len, max_len] = edgeStats( points, edgeG )
% [ avg_len, min_len, max_len] = edgeStats( points, edgeG )
%  statistics on edge lengths of vertex graph
%  points is Nx3 list of vertices
%  edgeG is (sparse) edge connectivity graph  (non-zero -> edge)

% [ RMS TODO AUGH ]
%   - this is wrong! interior edge lengths are included twice!

[i,j] = find( edgeG ~= 0 );
dists = sqrt( sum( (points(i,:)-points(j,:)).^2, 2 ) );

avg_len = mean(dists);
min_len = min(dists);
max_len = max(dists);


end
