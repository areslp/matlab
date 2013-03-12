function varargout = load_data_1_var(varargin)
% LOAD_DATA_1_VAR M-file for load_data_1_var.fig
%      LOAD_DATA_1_VAR, by itself, creates a new LOAD_DATA_1_VAR or raises the existing
%      singleton*.
%
%      H = LOAD_DATA_1_VAR returns the handle to a new LOAD_DATA_1_VAR or the handle to
%      the existing singleton*.
%
%      LOAD_DATA_1_VAR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LOAD_DATA_1_VAR.M with the given input arguments.
%
%      LOAD_DATA_1_VAR('Property','Value',...) creates a new LOAD_DATA_1_VAR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before load_data_1_var_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to load_data_1_var_OpeningFcn via varargin.
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

% Edit the above text to modify the response to help load_data_1_var

% Last Modified by GUIDE v2.5 23-Jul-2008 17:31:19

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @load_data_1_var_OpeningFcn, ...
                   'gui_OutputFcn',  @load_data_1_var_OutputFcn, ...
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


% --- Executes just before load_data_1_var is made visible.
function load_data_1_var_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to load_data_1_var (see VARARGIN)

% Choose default command line output for load_data_1_var
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes load_data_1_var wait for user response (see UIRESUME)
uiwait(handles.figure1);
% uiwait;



% --- Outputs from this function are returned to the command line.
function varargout = load_data_1_var_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function cn_Callback(hObject, eventdata, handles)
% hObject    handle to cn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cn as text
%        str2double(get(hObject,'String')) returns contents of cn as a double


% --- Executes during object creation, after setting all properties.
function cn_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cn (see GCBO)
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
if get(handles.i,'value')
    set(handles.cn,'Enable','on');
    set(handles.cnt,'ForegroundColor',[0 0 0]);
else
    set(handles.cn,'Enable','off');
    set(handles.cnt,'ForegroundColor',[0.4 0.4 0.4]);
end


