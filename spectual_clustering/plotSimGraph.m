function plotSimGraph( T, W, Colored, VertixCol, EdgeCol)
%PLOTSIMGRAPH Plots similarity graph
%   plotSimGraph(T, W) plots the similarity graph W of the data
%   set T, where T is a d-by-n matrix and W is the n-by-n graph.
%
%   'T' - Data points (d-by-n matrix)
%   'W' - Similarity graph (n-by-n matrix)
%   'Colored' - set this to one of the default colormaps (use 1
%      for the colormap 'Hot') to color the edges correspondingly
%      to their weight.
%   'VertixCol' (optional) - Color to be used to draw the
%      vertices (only use the default abbreviations like 'g'
%      etc.)
%   'EdgeCol' (optional) - In case Colored is not activated, 
%      define which color to use for plotting the edges (only use
%      the default abbrevations like 'g' etc.)

% catch illegal input
if nargin < 2
    error(['Not enough arguments: Data points and Similarity' ...
           ' Graph needed.']);
end

% if edges are to be colored, get input and set it up
if nargin > 2 && isnumeric(Colored) && Colored == 1
    Colored = 'Hot';
elseif nargin == 2
    Colored = 0;
end

% set up edge and vertix colors if not given
if nargin < 5
    EdgeCol = 'b';
    
    if nargin < 4
        VertixCol = 'r';
    end
end

% preallocate memory
xVal = zeros(2 * size(T, 2), 1);
yVal = zeros(2 * size(T, 2), 1);

% for colored edges we will have to do some calculations, so we
% set up some variables and prepare the colormap
if Colored ~= 0
    sMin = min(min(nonzeros(W)));
    sMax = max(max(W));
    if sMax == sMin
        sMin = 0;
    end
    cmap = colormap(Colored);
end

hold on;

% loop through all data points
for ii = 1:size(T, 2)
    
    % find vertices connected with the current vertix
    nnzs = find(W(ii, :) > 0);
    nnzs_ind = length(nnzs);
    
    if nnzs_ind > 0
        
        % if this vertix is connected to any else, store the
        % edges in the variables
        xVal(1 : 2 : 2 * nnzs_ind) = T(1, ii);
        xVal(2 : 2 : 2 * nnzs_ind) = T(1, nnzs);
        yVal(1 : 2 : 2 * nnzs_ind) = T(2, ii);
        yVal(2 : 2 : 2 * nnzs_ind) = T(2, nnzs);
        
        % plot edges
        if Colored == 0
            plot(xVal(1:2*nnzs_ind), yVal(1:2*nnzs_ind), ['-' ...
                 EdgeCol]);
        else
            for jj = 1:2:2*nnzs_ind
                curColVal = 1 - (W(ii, nnzs((jj+1)/2)) - ...
                                 sMin) / (sMax - sMin);
                curCol = cmap(1 + floor((size(cmap, 1)-1) * ...
                              curColVal), :);
                          
            plot(xVal(jj:jj+1), yVal(jj:jj+1), '-', ...
                 'Color', curCol);
            end
        end
    end
    
end

% plot vertices 
scatter(T(1, :), T(2, :), 5, VertixCol, 'filled');

hold off;

end

