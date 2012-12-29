%---------------------------
% Bachelor Thesis
% Author : Ingo Bürk
% Year   : 2011/2012
%---------------------------
% CLUSTERING OPTIONS
%---------------------------

% Choose similarity graph type
% - 2 for epsilon neighborhood
% - 3 for k-Nearest neighbors
% - 4 for mutual k-Nearest neighbors
SimGraphType   = 4;

% Choose spectral clustering algorithm
% - 1 for unnormalized Clustering
% - 2 for normalized   Clustering (Shi and Malik)
% - 3 for normalized   Clustering (Jordan and Weiss)
ClusteringType = 2;

% Choose a test mode
% - 1 for three Gaussians
% - 2 for two half moons
test_mode = 2;

% Choose whether to draw the similarity graph and if so,
% whether to color it (Only use with few data points!)
bplotSimGraph   = 0;
plotSimGraphMap = 1;

%---------------------------
% CREATE DATA
%---------------------------
time_data = tic;

if test_mode == 1
    
    % Three Gaussians
    T = [GenData_Gaussian(1000, [0 0], 0.1*eye(2)) ...
         GenData_Gaussian(1000, [2 0], 0.1*eye(2)) ...
         GenData_Gaussian(2000, [1 2], 0.2*eye(2))];
    
    num_clust = 3;
    Paramk    = 100;
    ParamEps  = 0.4;
    
elseif test_mode == 2
    
    % Two half moons
    T = [GenData_Ellipse(100, 1, 0.01, [0 0], 0, pi) ...
         GenData_Ellipse(100, 1, 0.01, [1 0.5], pi, 2*pi)];
    
    num_clust = 2;
    Paramk    = 15;
    ParamEps  = 0.4;
    
end

time_data = toc(time_data);
fprintf('Generating Data: %.2fs\n', time_data);

%---------------------------
% EXECUTE ALGORITHM AND PLOT
%---------------------------
 
% Parameter corresponding to SimGraphType (k or Eps)
switch SimGraphType
    case 2
        Param = ParamEps;
    case {3, 4}
        Param = Paramk;
end

% calculate adjacency matrix
time_graph = tic;

W = SimGraph(T, SimGraphType, Param, 1);

time_graph = toc(time_graph);
fprintf('Distance Matrix and Similarity Graph: %.2fs\n', ...
    time_graph);

% check for connected components (needs Bioinformatics Toolbox)
try
    conncomps = graphconncomp(W, 'Directed', false);
    if conncomps > 2
        warning(['Found %d connected components, maybe try ' ...
                 'changing the parameters.\n'], conncomps);
    else
        fprintf('Number of connected components: %d\n', ...
            conncomps);
    end
catch exception
    disp('Could not check for number of connected components.');
end

% compute clusters using spectral clustering
time_algorithm = tic;

A = SpectralClustering(W, num_clust, num_clust, ClusteringType);

time_algorithm = toc(time_algorithm);
fprintf('Clustering Algorithm: %.2fs\n', time_algorithm);

% plot result
scrsz = get(0,'ScreenSize');
hFig  = figure('Position', [scrsz(3)/10 scrsz(4)/3 ...
                            scrsz(3)/1.25 scrsz(4)/2]);
set(gca, 'Color', 'none');

if bplotSimGraph == 1
    subplot(1, 2, 1);
end
hold on;
col = [A zeros(size(A, 1), 1)];
sFigureTitle = sprintf(...
    'Clustered Data (Algorithm Time: %fs)', time_algorithm);
title(sFigureTitle)
cols = 'brgcmykw';

for ii = 1:num_clust
    curcol = [cols(mod(ii, size(cols, 2))) '*'];
    scatter(T(1, A(:, ii) == 1), T(2, A(:, ii) == 1), 7, curcol);
end
hold off;

% plot similarity graph
if bplotSimGraph == 1
    hSG = subplot(1, 2, 2);
    title('Similarity Graph');
    plotSimGraph(T, W, plotSimGraphMap);
end

% Commands to export with export_fig
%export_fig 'output.pdf' -q101 -transparent
%export_fig(hSG, 'output.pdf', '-q101', '-transparent')
