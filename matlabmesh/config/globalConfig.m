function [ value ] = globalConfig( name )
% [ value ] = globalConfig( setting )
%   returns global configuration value for string name
%     'enableMeshCaching' : 1/0 enable/disable mesh caching
%     'meshVersion' : mesh format version number  

if strcmp(name, 'meshVersion')
    value = 1;
end

switch name
    case 'meshVersion'
        value = 1.0;
    case 'pointsetVersion'    
        value = 1.1;

    case 'cachePath'
        thispath = which('globalConfig.m');
        [mat,idx] = regexp(thispath, '\\', 'match','end');
        cachepath = [ thispath(1:idx(numel(idx))), 'cache\' ];
        value = cachepath;
        
    case 'enableMeshCaching'
        value = 1;

    case 'enableWeightMatrixCaching'
        value = 0;
        
        
    otherwise
        error(['[globalConfig] unknown name ',name]);
end    


end
