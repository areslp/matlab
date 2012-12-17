function [hks] = HKS(filename, scale, opt, t)

% commond [hks] = HKS(filename, scale, opt)
% INPUTS
%  filename:  off file of triangle mesh.
%  opt.htype: the way to compute the parameter h. h = hs * neighborhoodsize
%             if htype = 'ddr' (data driven); h = hs if hytpe = 'psp' (pre-specify)
%             Default : 'ddr'
%  opt.hs:    the scaling factor that scales the neighborhood size to the
%             parameter h	where h^2 = 4t.
%             Default: 2, must > 0
%  opt.rho:   The cut-off for Gaussion function evaluation. 
%             Default: 3, must > 0
%  opt.dtype: the way to compute the distance 
%             dtype = 'euclidean' or 'geodesic';
%             Default : 'euclidean'
%  scale:  if scale = true, output the scaled hks
%          o.w. ouput the hks that is not scaled

% OUTPUTS
%  hks: ith row in this matrix is the heat kernel signature of the ith vertex


if nargin < 1
    error('Too few input arguments');	 
elseif nargin < 2
   scale = true
   opt.hs = 2;
	opt.rho = 3;
	opt.htype = 'ddr';
	opt.dtype = 'euclidean'

elseif nargin < 3
	opt.hs = 2;
	opt.rho = 3;
	opt.htype = 'ddr';
	opt.dtype = 'euclidean'
end
opt=parse_opt(opt)

if opt.hs <= 0 | opt.rho <= 0
	error('Invalid values in opt');
end

[W A h] = symmshlp_matrix(filename, opt);
Am = sparse([1:length(A)], [1:length(A)], A);

nev = min(300, length(A));

[evecs evals] = eigs(W, Am, nev, -1e-5);
evals = diag(evals)

%area = sum(A);
%A = (1/area) * A;
%evals = area * evals;
%evecs = sqrt(area) * evecs;
                     
tmin = abs(4*log(10) / evals(end));
tmax = abs(4*log(10) / evals(2));
nstep = 100;

stepsize = (log(tmax) - log(tmin)) / nstep;
logts = log(tmin):stepsize:log(tmax);
ts = exp(logts);
fprintf(1,'default ts is %f\n',ts);

if t~=0
	ts=t;
end

if scale == true, 
   hks = abs( evecs(:, 2:end) ).^2 * exp( ( abs(evals(2)) - abs(evals(2:end)) )  * ts);
   %Am = sparse([1:length(A)], [1:length(A)], A);
   colsum = sum(Am*hks);
   scale = 1.0./ colsum; 
   scalem = sparse([1:length(scale)], [1:length(scale)], scale);
   hks = hks * scalem;
else
   hks = abs( evecs(:, 2:end) ).^2 * exp( - abs(evals(2:end)) * ts);
end



% Parsing Option.
function option = parse_opt(opt)
option = opt;
option_names = {'hs', 'rho', 'htype', 'dtype'};
if ~isfield(option,'hs'),
	option = setfield(option,'hs',2);
end
if ~isfield(option,'rho'),
	option = setfield(option,'rho', 3);
end

if ~isfield(option,'htype'),
	option = setfield(option,'htype', 'ddr');
end

if ~isfield(option,'dtype'),
	option = setfield(option,'dtype', 'euclidean');
end


