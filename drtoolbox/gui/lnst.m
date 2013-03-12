% This file is part of the Matlab Toolbox for Dimensionality Reduction v0.7.2b.
% The toolbox can be obtained from http://homepage.tudelft.nl/19j49
% You are free to use, change, or redistribute this code in any way you
% want for non-commercial purposes. However, it is appreciated if you 
% maintain the name of the original author.
%
% (C) Laurens van der Maaten, 2010
% University California, San Diego / Delft University of Technology

function lnst(hs,he)
ls=get(he,'value');

switch ls
    case 1
        set(hs,'LineStyle','-','marker','x','color','b');
    case 2
        set(hs,'LineStyle','-','marker','o','color','r');
    case 3
        set(hs,'LineStyle','none','marker','o','color','b');
    case 4
        set(hs,'LineStyle','-','marker','none','color','g');
    case 5
        set(hs,'LineStyle','--','marker','d','color','b');
        
end
drawnow;