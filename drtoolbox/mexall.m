function mexall
%MEXALL Compiles all MEX-files of the Matlab Toolbox for Dimensionality Reduction
%
%   mexall
%
% Compiles all MEX-files of the Matlab Toolbox for Dimensionality Reduction.
%
%

% This file is part of the Matlab Toolbox for Dimensionality Reduction.
% The toolbox can be obtained from http://homepage.tudelft.nl/19j49
% You are free to use, change, or redistribute this code in any way you
% want for non-commercial purposes. However, it is appreciated if you 
% maintain the name of the original author.
%
% (C) Laurens van der Maaten, Delft University of Technology


    disp('Compiling...');
    cd techniques    
    try 
        mex -O computegr.c
    catch
        warning('Compiling failed. FastMVU might not work properly.');
    end
    try 
        mex -O mexCCACollectData2.c
    catch
        warning('Compiling failed. FastMVU might not work properly.');
    end
    try 
        mex -O mexCCACollectData.c
    catch
        warning('Compiling failed. CCA might not work properly.');
    end
    try
        if any(strcmpi(computer, {'MACI64', 'PCWIN64', 'GLNXA64', 'SOL64'}))
            mex -O -largeArrayDims dijkstra.cpp
        else
            mex -O dijkstra.cpp
        end            
    catch
        warning('Compiling failed. Isomap and LandmarkIsomap might not work properly.');
        if ispc
            warning('This error may be resolved by using the MinGW compiler instead of your current compiler.');
        end
    end
    cd ..
    disp('Compilation completed.');
    