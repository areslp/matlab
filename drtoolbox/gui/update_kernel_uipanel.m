% This file is part of the Matlab Toolbox for Dimensionality Reduction v0.7.2b.
% The toolbox can be obtained from http://homepage.tudelft.nl/19j49
% You are free to use, change, or redistribute this code in any way you
% want for non-commercial purposes. However, it is appreciated if you 
% maintain the name of the original author.
%
% (C) Laurens van der Maaten, 2010
% University California, San Diego / Delft University of Technology

if strcmpi(get(handles.uipanel_k,'Visible'),'on')
  
  if get(handles.kg,'value')
      set(handles.kgs,'visible','on');
      set(handles.kgst,'visible','on');
  else
      set(handles.kgs,'visible','off');
      set(handles.kgst,'visible','off');
  end
  
  if get(handles.kp,'value')
      
      set(handles.kpR,'visible','on');
      set(handles.kpRt,'visible','on');
      
      set(handles.kpd,'visible','on');
      set(handles.kpdt,'visible','on');
      
  else
      
      set(handles.kpR,'visible','off');
      set(handles.kpRt,'visible','off');
      
      set(handles.kpd,'visible','off');
      set(handles.kpdt,'visible','off');
      
  end
else
    
  set(handles.kgs,'visible','off');
  set(handles.kgst,'visible','off');
  
  set(handles.kpR,'visible','off');
  set(handles.kpRt,'visible','off');
  
  set(handles.kpd,'visible','off');
  set(handles.kpdt,'visible','off');
  
end
      
    