function varargout = choose_method(varargin)
% CHOOSE_METHOD M-file for choose_method.fig
%      CHOOSE_METHOD, by itself, creates a new CHOOSE_METHOD or raises the existing
%      singleton*.
%
%      H = CHOOSE_METHOD returns the handle to a new CHOOSE_METHOD or the handle to
%      the existing singleton*.
%
%      CHOOSE_METHOD('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CHOOSE_METHOD.M with the given input arguments.
%
%      CHOOSE_METHOD('Property','Value',...) creates a new CHOOSE_METHOD or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before choose_method_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to choose_method_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% This file is part of the Matlab Toolbox for Dimensionality Reduction v0.7.2b.
% The toolbox can be obtained from http://homepage.tudelft.nl/19j49
% You are free to use, change, or redistribute this code in any way you
% want for non-commercial purposes. However, it is appreciated if you 
% maintain the name of the original author.
%
% (C) Laurens van der Maaten, 2010
% University California, San Diego / Delft University of Technology

% Edit the above text to modify the response to help choose_method

% Last Modified by GUIDE v2.5 25-Jul-2008 10:40:42

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @choose_method_OpeningFcn, ...
                   'gui_OutputFcn',  @choose_method_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before choose_method is made visible.
function choose_method_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to choose_method (see VARARGIN)

% Choose default command line output for choose_method
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes choose_method wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = choose_method_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
vl=get(hObject,'Value');
switch vl
    case 1
        str='CorrDim';
    case 2
        str='NearNbDim';
    case 3
        str='GMST';
    case 4
        str='PackingNumbers';
    case 5
        str='EigValue';
    case 6
        str='MLE';
end
set(handles.str,'string',str);
drawnow;
        


% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function str_Callback(hObject, eventdata, handles)
% hObject    handle to str (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of str as text
%        str2double(get(hObject,'String')) returns contents of str as a double


% --- Executes during object creation, after setting all properties.
function str_CreateFcn(hObject, eventdata, handles)
% hObject    handle to str (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ok.
function ok_Callback(hObject, eventdata, handles)
uiresume(handles.figure1);


