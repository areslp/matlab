function varargout = mapping_parameters(varargin)
% MAPPING_PARAMETERS M-file for mapping_parameters.fig
%      MAPPING_PARAMETERS, by itself, creates a new MAPPING_PARAMETERS or raises the existing
%      singleton*.
%
%      H = MAPPING_PARAMETERS returns the handle to a new MAPPING_PARAMETERS or the handle to
%      the existing singleton*.
%
%      MAPPING_PARAMETERS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MAPPING_PARAMETERS.M with the given input arguments.
%
%      MAPPING_PARAMETERS('Property','Value',...) creates a new MAPPING_PARAMETERS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before mapping_parameters_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to mapping_parameters_OpeningFcn via varargin.
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

% Last Modified by GUIDE v2.5 28-Jul-2008 16:23:20

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @mapping_parameters_OpeningFcn, ...
                   'gui_OutputFcn',  @mapping_parameters_OutputFcn, ...
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


% --- Executes just before mapping_parameters is made visible.
function mapping_parameters_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to mapping_parameters (see VARARGIN)

% Choose default command line output for mapping_parameters
handles.output = hObject;

no_dims=varargin{1};
handles.islb=varargin{2};

if length(no_dims)~=0
    set(handles.nd,'string',num2str(round(no_dims)));
else
    set(handles.nd,'string','2');
end

case1; % case one is default for start

% Update handles structure
guidata(hObject, handles);

% handles

% UIWAIT makes mapping_parameters wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = mapping_parameters_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
%     figure1: 218.0023
%             sm: 240.0021
%             na: 239.0021
%            nat: 238.0021
%     uipanel_tp: 235.0021
%            prp: 234.0021
%           prpt: 233.0021
%            kpd: 232.0021
%           kpdt: 231.0021
%            kpR: 230.0021
%           kpRt: 229.0021
%            kgs: 228.0021
%           kgst: 227.0021
%      uipanel_k: 223.0022
%     uipanel_ft: 220.0022
%            sig: 53.0027
%           sigt: 52.0027
%     uipanel_ei: 49.0027
%            prc: 48.0027
%           prct: 47.0027
%              k: 46.0028
%             kt: 45.0028
%             mi: 44.0028
%            mit: 43.0029
%             nd: 42.0029
%          text3: 41.0034
%             wl: 40.0029
%          text1: 39.0033
%       listbox1: 219.0023
%             tl: 237.0021
%             tg: 236.0021
%             kp: 226.0022
%             kg: 225.0022
%             kl: 224.0022
%            ftn: 222.0022
%            fty: 221.0022
%            eij: 51.0027
%            eim: 50.0027
%         output: 218.0023
%           islb:

% PCA
% LDA
% MDS
% SimplePCA
% ProbPCA
% FactorAnalysis
% Isomap
% LandmarkIsomap
% LLE
% Laplacian
% HessianLLE
% LTSA
% MVU
% CCA
% LandmarkMVU
% FastMVU
% DiffusionMaps
% KernelPCA
% GDA
% SNE
% SymSNE
% t-SNE
% LPP
% NPE
% LLTSA
% SPE
% Autoencoder
% LLC
% ManifoldChart
% CFA
% GPLVM

switch get(hObject,'Value')
    case 1 % PCA
        % no parameters;
         case1;
    case 2 % LDA
        % no parameters;
         case1;
         if handles.islb
             set(handles.wl,'visible','off');
             set(handles.sm,'Enable','on');
         else
             set(handles.wl,'visible','on');
             set(handles.sm,'Enable','off');
         end
    case 3 % MDS
        % no parameters;
        case1;
    case 4 % SimplePCA
        % no parameters;
        case1;
    case 5 % ProbPCA
        case1; % hide all contols first
        set(handles.mi,'Visible','on');
        set(handles.mit,'Visible','on');
    case 6 % FactorAnalysis
        % no parameters;
        case1;
    case 7 % Isomap
        case1; % hide all contols first
        set(handles.k,'Visible','on','string',12);
        set(handles.kt,'Visible','on');
        set(handles.ka,'Visible','on');
        adaptive_callback;
    case 8 % LandmarkIsomap
        case1; % hide all contols first
        set(handles.k,'Visible','on','string',12);
        set(handles.kt,'Visible','on');
        set(handles.ka,'Visible','on');
        adaptive_callback;
        
        set(handles.prc,'Visible','on');
        set(handles.prct,'Visible','on');
    case 9 % LLE
        case1; % hide all contols first
        set(handles.k,'Visible','on','string',12);
        set(handles.kt,'Visible','on');
        set(handles.ka,'Visible','on');
        adaptive_callback;
        
        set(handles.uipanel_ei,'Visible','on');
    case 10 % Laplacian
        case1; % hide all contols first
        set(handles.k,'Visible','on','string',12);
        set(handles.kt,'Visible','on');
        set(handles.ka,'Visible','on');
        adaptive_callback;
        
        set(handles.sig,'Visible','on');
        set(handles.sigt,'Visible','on');
        
        set(handles.uipanel_ei,'Visible','on');
    case 11 % HessianLLE
        case1; % hide all contols first
        set(handles.k,'Visible','on','string',12);
        set(handles.kt,'Visible','on');
        set(handles.ka,'Visible','on');
        adaptive_callback;
        
        set(handles.uipanel_ei,'Visible','on');
    case 12 % LTSA
        case1; % hide all contols first
        set(handles.k,'Visible','on','string',12);
        set(handles.kt,'Visible','on');
        set(handles.ka,'Visible','on');
        adaptive_callback;
        
        set(handles.uipanel_ei,'Visible','on');
    case 13 % MVU
        case1; % hide all contols first
        set(handles.k,'Visible','on','string',12);
        set(handles.kt,'Visible','on');
        set(handles.ka,'Visible','on');
        adaptive_callback;
        
        set(handles.uipanel_ei,'Visible','on');
    case 14 % CCA
        case1; % hide all contols first
        set(handles.k,'Visible','on','string',12);
        set(handles.kt,'Visible','on');
        set(handles.ka,'Visible','on');
        adaptive_callback;
        
        set(handles.uipanel_ei,'Visible','on');
    case 15 % LandmarkMVU
        case1; % hide all contols first
        set(handles.k,'Visible','on','string',5);
        set(handles.kt,'Visible','on');
        set(handles.ka,'Visible','on');
        adaptive_callback;
    case 16 % FastMVU
        case1; % hide all contols first
        set(handles.k,'Visible','on','string',5);
        set(handles.kt,'Visible','on');
        set(handles.ka,'Visible','on');
        adaptive_callback;
        
        set(handles.uipanel_ei,'Visible','on');
        
        set(handles.uipanel_ft,'Visible','on');
    case 17 % DiffusionMaps
        case1; % hide all contols first
        set(handles.t,'Visible','on');
        set(handles.tt,'Visible','on');
        
        set(handles.sig,'Visible','on');
        set(handles.sigt,'Visible','on');
    case 18 % KernelPCA
        case1; % hide all contols first
                
        set(handles.uipanel_k,'Visible','on');
        update_kernel_uipanel;
    case 19 % GDA
        case1; % hide all contols first
        if handles.islb
             set(handles.wl,'visible','off');
             set(handles.sm,'Enable','on');
             
             set(handles.uipanel_k,'Visible','on');
             update_kernel_uipanel;
             
         else
             set(handles.wl,'visible','on');
             set(handles.sm,'Enable','off');
        end
    case 20 % SNE
        case1; % hide all contols first
        set(handles.prp,'Visible','on');
        set(handles.prpt,'Visible','on');
    case 21 % SymSNE
        case1; % hide all contols first
        set(handles.prp,'Visible','on');
        set(handles.prpt,'Visible','on');
    case 22 % t-SNE
        case1; % hide all contols first
        set(handles.prp,'Visible','on');
        set(handles.prpt,'Visible','on');
    case 23 % LPP
        case1; % hide all contols first
        set(handles.k,'Visible','on','string',12);
        set(handles.kt,'Visible','on');
        set(handles.ka,'Visible','on');
        adaptive_callback;
        
        set(handles.sig,'Visible','on');
        set(handles.sigt,'Visible','on');
        
        set(handles.uipanel_ei,'Visible','on');
    case 24 % NPE
        case1; % hide all contols first
        set(handles.k,'Visible','on','string',12);
        set(handles.kt,'Visible','on');
        set(handles.ka,'Visible','on');
        adaptive_callback;
       
        set(handles.uipanel_ei,'Visible','on');
    case 25 % LLTSA
        case1; % hide all contols first
        set(handles.k,'Visible','on','string',12);
        set(handles.kt,'Visible','on');
        set(handles.ka,'Visible','on');
        adaptive_callback;
       
        set(handles.uipanel_ei,'Visible','on');
    case 26 % SPE
        case1; % hide all contols first
        set(handles.uipanel_tp,'Visible','on');
        update_type_uipanel;
        set(handles.k,'Enable','on');
        adaptive_callback;
    case 27 % Autoencoder
        % no parameters;
         case1;
    case 28 % LLC
        case1; % hide all contols first
        set(handles.k,'Visible','on','string',12);
        set(handles.kt,'Visible','on');
        set(handles.ka,'Visible','on');
        adaptive_callback;
        
        set(handles.na,'Visible','on');
        set(handles.nat,'Visible','on');
        
        set(handles.mi,'Visible','on');
        set(handles.mit,'Visible','on');
        
        set(handles.uipanel_ei,'Visible','on');
    case 29 % ManifoldChart
        case1; % hide all contols first
        
        set(handles.na,'Visible','on');
        set(handles.nat,'Visible','on');
        
        set(handles.mi,'Visible','on');
        set(handles.mit,'Visible','on');
        
        set(handles.uipanel_ei,'Visible','on');
    case 30 % CFA
        case1; % hide all contols first
        
        set(handles.na,'Visible','on');
        set(handles.nat,'Visible','on');
        
        set(handles.mi,'Visible','on');
        set(handles.mit,'Visible','on');
    case 31 % GPLVM
        case1; % hide all contols first
        
        set(handles.sig,'Visible','on');
        set(handles.sigt,'Visible','on');
    case 32 % NCA
        case1;
        
    case 33 % MCML
        case1;
             
end


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



function nd_Callback(hObject, eventdata, handles)
% hObject    handle to nd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nd as text
%        str2double(get(hObject,'String')) returns contents of nd as a double


% --- Executes during object creation, after setting all properties.
function nd_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function mi_Callback(hObject, eventdata, handles)
% hObject    handle to mi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mi as text
%        str2double(get(hObject,'String')) returns contents of mi as a double


% --- Executes during object creation, after setting all properties.
function mi_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function k_Callback(hObject, eventdata, handles)
% hObject    handle to k (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of k as text
%        str2double(get(hObject,'String')) returns contents of k as a double


% --- Executes during object creation, after setting all properties.
function k_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function prc_Callback(hObject, eventdata, handles)
% hObject    handle to prc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of prc as text
%        str2double(get(hObject,'String')) returns contents of prc as a double


% --- Executes during object creation, after setting all properties.
function prc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to prc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sig_Callback(hObject, eventdata, handles)
% hObject    handle to sig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sig as text
%        str2double(get(hObject,'String')) returns contents of sig as a double


% --- Executes during object creation, after setting all properties.
function sig_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function kgs_Callback(hObject, eventdata, handles)
% hObject    handle to kgs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of kgs as text
%        str2double(get(hObject,'String')) returns contents of kgs as a double


% --- Executes during object creation, after setting all properties.
function kgs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to kgs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function kpR_Callback(hObject, eventdata, handles)
% hObject    handle to kpR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of kpR as text
%        str2double(get(hObject,'String')) returns contents of kpR as a double


% --- Executes during object creation, after setting all properties.
function kpR_CreateFcn(hObject, eventdata, handles)
% hObject    handle to kpR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function kpd_Callback(hObject, eventdata, handles)
% hObject    handle to kpd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of kpd as text
%        str2double(get(hObject,'String')) returns contents of kpd as a double


% --- Executes during object creation, after setting all properties.
function kpd_CreateFcn(hObject, eventdata, handles)
% hObject    handle to kpd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function prp_Callback(hObject, eventdata, handles)
% hObject    handle to prp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of prp as text
%        str2double(get(hObject,'String')) returns contents of prp as a double


% --- Executes during object creation, after setting all properties.
function prp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to prp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function na_Callback(hObject, eventdata, handles)
% hObject    handle to na (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of na as text
%        str2double(get(hObject,'String')) returns contents of na as a double


% --- Executes during object creation, after setting all properties.
function na_CreateFcn(hObject, eventdata, handles)
% hObject    handle to na (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in sm.
function sm_Callback(hObject, eventdata, handles)
set(handles.sm,'Enable','off');
drawnow;
uiresume(handles.figure1);




function t_Callback(hObject, eventdata, handles)
% hObject    handle to t (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of t as text
%        str2double(get(hObject,'String')) returns contents of t as a double


% --- Executes during object creation, after setting all properties.
function t_CreateFcn(hObject, eventdata, handles)
% hObject    handle to t (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --------------------------------------------------------------------
function uipanel_k_SelectionChangeFcn(hObject, eventdata, handles)
update_kernel_uipanel;




% --------------------------------------------------------------------
function uipanel_tp_SelectionChangeFcn(hObject, eventdata, handles)
update_type_uipanel;




% --- Executes on button press in ka.
function ka_Callback(hObject, eventdata, handles)
adaptive_callback;


