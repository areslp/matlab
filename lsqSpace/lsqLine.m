function [Point, Dir,varargout] = lsqLine(V,varargin)
    % [Point, Dir] = lsqLine(V)
    %
    % Given an n x m matrix V, whose rows are a set of sample m-vectors,
    % the routine calculates a line in R^m which optimally represents the samples in a least-square 
    % sense. The line is returned as a Point, and Dir row vectors.
    %
    % [...] = lsqLine(...,'discardworst',k)
    % Iteratively discards the k samples farthest from the optimal line found.
    %
    % [...] = lsqLine(...,'discardthresh',th)
    % Iteratively discards samples whose distance from the found line exceeds th .
    %
    % [... , samplesused] = lsqLine(...)
    % Returns a vector of indices of the samples that weren't discarded in the process.
    %

    
    nin = nargin - 1;
    nout = nargout - 2;

    varargout = {};
    
    if nout>1
	  error('Improper number of output arguments in lsqLine.');
    end

    if nin && nout
	  [Point, Dir ,varargout{1}] = lsqAffineSpace(V, 1,  varargin);
    elseif nin
	  [Point, Dir] = lsqAffineSpace(V, 1, varargin);
    else
	  [Point, Dir] = lsqAffineSpace(V, 1);
    end
