function [ Energy ] = FnOptLaplacianScale( paramtype, surface, k, W, varargin )
%FNOPTLAPLACIANSCALE Summary of this function goes here
%   Detailed explanation goes here

if strcmp(paramtype, 'lle')
    surface.u = embedLLE(surface.v, k*W, surface.vidx(surface.isboundaryv~=0), surface.loops, surface );
    [ EConformal, EDirichlet, Area, Turning ] = eConformal( surface, k*W );

elseif strcmp(paramtype, 'dncp')
    surface.u = embedDNCP(surface, k*W, varargin{:} );
    [ EConformal, EDirichlet, Area, Turning ] = eConformal( surface, k*W );
    
elseif strcmp(paramtype, 'scp')
    surface.u = embedSCP(surface, 'generalized', k*W );
    [ EConformal, EDirichlet, Area, Turning ] = eConformal( surface, k*W );
    
elseif strcmp(paramtype, 'pscp')
    surface.u = embedPSCP(surface, k*W);
    [ EConformal, EDirichlet, Area, Turning ] = eConformal( surface, k*W );
    
end

%scale = 1 / Area;
scale = 1;
Area = Area*scale;
EDirichlet = EDirichlet*scale;
EConformal = EDirichlet - Area;

%if EConformal < 0
%    EConformal = EConformal+10;
%end
    

% return this as energy...not continuous, but smooth between big jumps...
%Energy = EConformal + 10*(Turning-1);
%Energy = log(EConformal.^2) + 10*(Turning-1);
Energy = EConformal.^2 + 10*abs(Turning-1);
%Energy = abs(EConformal) + 10*(Turning-1);
   
fprintf('E: %3.10f  k: %3.5f C:%4.5f  =  D:%3.8f - A:%3.8f  (turn %2d)\n', Energy, k, EConformal, EDirichlet, Area, Turning);
    
end
