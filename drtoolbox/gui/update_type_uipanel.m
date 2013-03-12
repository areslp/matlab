% This file is part of the Matlab Toolbox for Dimensionality Reduction v0.7.2b.
% The toolbox can be obtained from http://homepage.tudelft.nl/19j49
% You are free to use, change, or redistribute this code in any way you
% want for non-commercial purposes. However, it is appreciated if you 
% maintain the name of the original author.
%
% (C) Laurens van der Maaten, 2010
% University California, San Diego / Delft University of Technology

if strcmpi(get(handles.uipanel_tp,'Visible'),'on')
    if get(handles.tl,'value')
        set(handles.k,'Visible','on','string',12);
        set(handles.kt,'Visible','on');
        %set(handles.ka,'Visible','on');
    else
        set(handles.k,'Visible','off');
        set(handles.kt,'Visible','off');
        %set(handles.ka,'Visible','off');
    end
end