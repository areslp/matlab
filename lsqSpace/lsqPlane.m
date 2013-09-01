function [Point, Normal,varargout] = lsqPlane(V,varargin)
    % [Point, Normal] = lsqPlane(V)
    %
    % Given an n x m matrix V, whose rows are a set of sample m-vectors,
    % the routine calculates a hyperplane in R^m which optimally represents the samples 
    % in a least-square sense. The hyperplane is returned as a Point and Normal row 
    % vectors.
    %
    % [...] = lsqPlane(...,'discardworst',k)
    % Iteratively discards the k samples farthest from the optimal plane found.
    %
    % [...] = lsqPlane(...,'discardthresh',th)
    % Iteratively discards samples whose distance from the found plane exceeds th .
    %
    % [... , samplesused] = lsqPlane(...)
    % Returns a vector of indices of the samples that weren't discarded in the process.
    %

    nin = nargin - 1;
    nout = nargout - 2;
    m=size(V,2);

    varargout = {};

    if nout>1
	  error('Improper number of output arguments in lsqPlane.');
    end

    if nin && nout
	  [Point, Normal ,varargout{1}] = lsqAffineSpace(V, m-1, 'orthogonal', varargin);
    elseif nin
	  [Point, Normal] = lsqAffineSpace(V, m-1, 'orthogonal', varargin);
    else
	  [Point, Normal] = lsqAffineSpace(V, m-1, 'orthogonal');
    end

