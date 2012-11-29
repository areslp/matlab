function [ P ] = genpoints( N, dim, type )
%[ P ] = genpoints( N, dim, type )
%   generate point samples
%     N = number of points 
%        !! approximate - may get fewer points depending on method
%     dim = dimension (some methods only work for d=2 or d=3
%     type:
%       'grid_unit_square' - regular-grid in [0,1]^(dim==2)
%       'stratified_unit_square' - jittered 'grid_unit_square'
%       'uniform_unit_square' - uniform random samples in [0,1]^dim

if dim == 2
    
    if strcmp(type, 'grid_unit_square')
        k = floor(sqrt(N));
        d = 1 / (k-1);
        x = repmat((0:k-1)*d,1,k);
        y = reshape( repmat((0:k-1)*d,k,1), k*k, 1);
        P = [x',y];
        
    elseif strcmp(type, 'stratified_unit_square')
        k = floor(sqrt(N));
        d = 1 / (k-1);
        x = repmat((0:k-1)*d,1,k);
        y = reshape( repmat((0:k-1)*d,k,1), k*k, 1);
        P = [x',y] + rand(k*k,2)*d/2;
        
    elseif strcmp(type, 'uniform_unit_square')
        P = rand(N,2);
    end
    
    
end


end
