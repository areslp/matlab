function recX = reconstruct_data(mappedX, mapping)
%RECONSTRUCT_DATA Reconstructs data from low-dimensional data representation
%
%   recX = reconstruct_data(mappedX, mapping)
%
% The function tries to reconstruct the original data from the
% low-dimensional representation in mappedX. The function requires the
% original mapping as input, too. This function only works for linear
% dimensionality reduction techniques and for autoencoders. (Other
% techniques do not support such backprojections.)
% The reconstructed data can be obtained from recX. You can use measure the
% SSE between the original X and recX to get a measure for the quality of
% the reconstruction.
%
%

% This file is part of the Matlab Toolbox for Dimensionality Reduction.
% The toolbox can be obtained from http://homepage.tudelft.nl/19j49
% You are free to use, change, or redistribute this code in any way you
% want for non-commercial purposes. However, it is appreciated if you 
% maintain the name of the original author.
%
% (C) Laurens van der Maaten, Delft University of Technology


    % Handle PRTools dataset
    if strcmp(class(mappedX), 'dataset')
        prtools = true;
        A = mappedX;
        mappedX = mappedX.data;        
    else 
        prtools = false;
    end

    % Perform reconstruction
    switch mapping.name
        
        case {'PCA', 'LDA', 'LPP', 'NPE', 'LLTSA', 'SPCA', 'PPCA', 'FA', 'NCA', 'MCML', 'LMNN'}
            recX = bsxfun(@plus, mappedX * mapping.M', mapping.mean);
            
        case {'Autoencoder'}
            recX = recon_data_from_autoenc(mapping.network, mappedX);
            
        otherwise
            error('Reconstruction of data is not supported for this technique.');
    end
    
    % Handle PRTools dataset
    if prtools
        A.data = recX;
        recX = A;
    end    
        