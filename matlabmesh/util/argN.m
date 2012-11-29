function [ E ] = argN( return_N, f, varargin )
%[ E ] = argN( return_N, f, varargin )
%   return Nth returned argument of f(varargin)

if return_N == 1
    E = f(varargin{:});
elseif return_N == 2
    [a,b] = f(varargin{:});
    cells = {a,b};
    E = cells{return_N};
else
    fprintf('argN does not support that number of arguments (add it yourself!)\n');
end
