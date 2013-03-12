% This file is part of the Matlab Toolbox for Dimensionality Reduction v0.7.2b.
% The toolbox can be obtained from http://homepage.tudelft.nl/19j49
% You are free to use, change, or redistribute this code in any way you
% want for non-commercial purposes. However, it is appreciated if you 
% maintain the name of the original author.
%
% (C) Laurens van der Maaten, 2010
% University California, San Diego / Delft University of Technology

if strcmpi(get(handles.k,'Visible'),'on')
    if get(handles.listbox1,'Value')~=25
        if get(handles.ka,'Value')
            set(handles.k,'Enable','off');
        else
            set(handles.k,'Enable','on');
        end
    end
end