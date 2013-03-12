function [x,obj] = hillclimber2c(DD, x, varargin)
%HILLCLIMBER Performs hillclimbing using initial solution 
%
%   function [x, obj] = hillclimber(DD, x, varargin)
%
%
%

% This file is part of the Matlab Toolbox for Dimensionality Reduction.
% The toolbox can be obtained from http://homepage.tudelft.nl/19j49
% You are free to use, change, or redistribute this code in any way you
% want for non-commercial purposes. However, it is appreciated if you 
% maintain the name of the original author.
%
% (C) Laurens van der Maaten, Delft University of Technology


    n = length(DD);
    
    % Set parameters for the hillclimber
    pars.maxiter = 10000;
    pars.stepsize = 1e-05;
    pars.ETA = 1e0-5;
    pars.truex = [];
    pars.verbose = 1;
    pars.acc = 1.01;
    pars.printevery = 100;
    pars.othresh = 1e-10;
    pars.eta = 0;

    % Initialize variables
    ii = find(DD);
    dd = full(DD(ii));
    [i2,i1] = ind2sub(size(DD), ii);
    dims = size(x);

    % Perform hill-climbing
    [x, obj, i] = minimize(x(:), 'hill_obj', -pars.maxiter, dims, [i1 i2], dd, pars);
    x = reshape(x, dims);