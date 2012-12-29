function [ T ] = GenData_Ellipse( num, r, var, pos, ...
    start_angle, end_angle )
%GENDATA_ELLIPSE Calculates a ellipse sample data set
%   GenData_Ellipse returns a matrix containing data points as
%   columns, which are distributed in the shape of an ellipse
%   with specified radii r and center point pos. Moreover, the
%   data points will be spread with variance var and start_angle
%   and end_angle define begin and start of what part of the
%   ellipse to use.
%
%   For example, GenData_Ellipse(1000, 1, 0, [0 0], 0, 2*pi) will
%   give 1000 data points in perfect shape of the unit circle.
%
%   'num' - Number of data points to be created
%   'r' - Either vector containing the radii or a scalar radius
%   'var' - scalar variance value to spread data
%   'pos' - vector defining the center of the ellipse
%   'start_angle' - start angle of the ellipse (radian)
%   'end_angle' - end angle of the ellipse (radian)
%
%   Author: Ingo Buerk
%   Year  : 2011/2012
%   Bachelor Thesis

% make column vectors
pos = pos(:);
r   = r(:);

% if scalar radius was given, make vector
if size(r, 1) == 1
    r = [r; r];
end

% create vector of angles
angles = rand(num, 1).*(end_angle-start_angle)+start_angle;

% create data set
T = zeros(2, num);
for i = 1:num
    T(1, i) = pos(1) + r(1) * cos(angles(i));
    T(2, i) = pos(2) + r(2) * sin(angles(i));
end

% spread data
T = T + sqrt(var) * randn(2, num);

end

