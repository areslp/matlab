function varargout = load_data_vars(varargin)
% LOAD_DATA_VARS M-file for load_data_vars.fig
%      LOAD_DATA_VARS, by itself, creates a new LOAD_DATA_VARS or raises the existing
%      singleton*.
%
%      H = LOAD_DATA_VARS returns the handle to a new LOAD_DATA_VARS or the handle to
%      the existing singleton*.
%
%      LOAD_DATA_VARS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LOAD_DATA_VARS.M with the given input arguments.
%
%      LOAD_DATA_VARS('Property','Value',...) creates a new LOAD_DATA_VARS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before load_data_vars_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to load_data_vars_OpeningFcn via varargin.
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

% Edit the above text to modify the response to help load_data_vars

% Last Modified by GUIDE v2.5 23-Jul-2008 18:19:22

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @load_data_vars_OpeningFcn, ...
                   'gui_OutputFcn',  @load_data_vars_OutputFcn, ...
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


% --- Executes just before load_data_vars is made visible.
function load_data_vars_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to load_data_vars (see VARARGIN)

% Choose default command line output for load_data_vars
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);


vn=varargin{1}; % variable names

str='';
for vnc=1:length(vn)
    str=[str vn{vnc}];
    if vnc~=length(vn)
        str=[str '   '];
    end
end

set(handles.lv,'string',str);

set(handles.d,'string',vn{1});
set(handles.ivn,'string',vn{2});


% UIWAIT makes load_data_vars wait for user response (see UIRESUME)
uiwait(handles.figure1);
%uiwait;


% --- Outputs from this function are returned to the command line.
function varargout = load_data_vars_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function lv_Callback(hObject, eventdata, handles)
% hObject    handle to lv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lv as text
%        str2double(get(hObject,'String')) returns contents of lv as a double


% --- Executes during object creation, after setting all properties.
function lv_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function d_Callback(hObject, eventdata, handles)
% hObject    handle to d (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of d as text
%        str2double(get(hObject,'String')) returns contents of d as a double


% --- Executes during object creation, after setting all properties.
function d_CreateFcn(hObject, eventdata, handles)
% hObject    handle to d (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ivn_Callback(hObject, eventdata, handles)
% hObject    handle to ivn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ivn as text
%        str2double(get(hObject,'String')) returns contents of ivn as a double


% --- Executes during object creation, after setting all properties.
function ivn_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ivn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function icn_Callback(hObject, eventdata, handles)
% hObject    handle to icn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of icn as text
%        str2double(get(hObject,'String')) returns contents of icn as a double


% --- Executes during object creation, after setting all properties.
function icn_CreateFcn(hObject, eventdata, handles)
% hObject    handle to icn (see GCBO)
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


% --------------------------------------------------------------------
function uipanel1_SelectionChangeFcn(hObject, eventdata, handles)

if get(handles.iv,'value')
    set(handles.ivn,'Enable','on');
else
    set(handles.ivn,'Enable','off');
end

if get(handles.ic,'value')
    set(handles.icn,'Enable','on');
else
    set(handles.icn,'Enable','off');
end



