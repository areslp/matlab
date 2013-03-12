function welcome
%WELCOME Displays DR Toolbox version information
%
%   welcome
%
% Displays DR Toolbox version information.

% This file is part of the Matlab Toolbox for Dimensionality Reduction.
% The toolbox can be obtained from http://homepage.tudelft.nl/19j49
% You are free to use, change, or redistribute this code in any way you
% want for non-commercial purposes. However, it is appreciated if you 
% maintain the name of the original author.
%
% (C) Laurens van der Maaten, Delft University of Technology

    global DR_WELCOME;
    if isempty(DR_WELCOME)
        disp(' ');
        disp('   Welcome to the Matlab Toolbox for Dimensionality Reduction, version 0.8b (18-April-2012).');
        disp('      You are free to modify or redistribute this code (for non-commercial purposes), as long as a reference');
        disp('      to the original author (Laurens van der Maaten, Delft University of Technology) is retained.');
        disp('      For more information, please visit http://homepage.tudelft.nl/19j49');
        disp(' ');
        DR_WELCOME = 1;
    end
