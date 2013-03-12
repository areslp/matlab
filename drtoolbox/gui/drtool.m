function varargout = drtool(varargin)
% DRTOOL M-file for drtool.fig
%      DRTOOL, by itself, creates a new DRTOOL or raises the existing
%      singleton*.
%
%      H = DRTOOL returns the handle to a new DRTOOL or the handle to
%      the existing singleton*.
%
%      DRTOOL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DRTOOL.M with the given input arguments.
%
%      DRTOOL('Property','Value',...) creates a new DRTOOL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before drtool_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to drtool_OpeningFcn via varargin.
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

% Edit the above text to modify the response to help drtool

% Last Modified by GUIDE v2.5 16-Sep-2008 12:18:29

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @drtool_OpeningFcn, ...
                   'gui_OutputFcn',  @drtool_OutputFcn, ...
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


% --- Executes just before drtool is made visible.
function drtool_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to drtool (see VARARGIN)

% Choose default command line output for drtool
handles.output = hObject;

set(handles.edr,'UserData',[]); % demention not estimated at begining

% here data will be stored:
handles.X=[];
handles.labels=[];
handles.islb=false; % if labels provided
handles.isl=false; % if data loaded
handles.mcd=false; % if mapping was calculated
handles.mX=[]; % mapped X

handles.m1={}; % m-code from part1
handles.m2={}; % m-code from part2
handles.m21={}; % addition with no_dims
handles.m3={}; % m-code from part3
handles.a2=false; % if code part with no_dims was writed
handles.mstf={}; % save to file part

handles.isxls=[]; % if data loaded as xls-file (will be used before save history)

%handles.ndr=[]; % rounded number of dimention from intrinsic_dim, need for m-file

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes drtool wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = drtool_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in ld.
function ld_Callback(hObject, eventdata, handles)
handles.a2=false;
set(handles.ld,'Enable','off');
drawnow;
% load data lead to clear mapped data
handles.mcd=false; % if mapping was calculated
handles.mX=[]; % mapped X
set(handles.cd,'Visible','off');

try
    hs=load_data;
catch
    set(handles.ld,'Enable','on');
    return
end
hands_l_d=guidata(hs);
if (~get(hands_l_d.file,'value'))&&(~get(hands_l_d.xls,'value'))
    %if not from file and not from xls-file
    handles.isxls=[];
    
    set(hands_l_d.ok,'Enable','off');
    drawnow;
    
    if get(hands_l_d.sw,'value')
        dataname='swiss';
    end
    
    if get(hands_l_d.hl,'value')
        dataname='helix';
    end
    
    if get(hands_l_d.tp,'value')
        dataname='twinpeaks';
    end
    
    if get(hands_l_d.cl,'value')
        dataname='3d_clusters';
    end
    
    if get(hands_l_d.is,'value')
        dataname='intersect';
    end
    
    n=str2num(get(hands_l_d.ndp,'string'));
    noise=str2num(get(hands_l_d.nl,'string'));
    
    [X, labels] = generate_data(dataname, n,noise);
    handles.X=X;
    handles.labels=labels;
    handles.islb=true;
    handles.isl=true;
    set(handles.dld,'visible','on');
    
    % m-file memorizing:
    % clear all previose history because new original dataset generated:
    handles.m1={}; 
    handles.m2={}; 
    handles.m3={};
    handles.mstf={};
    
    handles.m1={ '% generate data';
                ['n = ' num2str(n) '; % number of datapoints' ];
                ['noise = ' num2str(noise) '; % noise level' ];
                ['[X, labels] = generate_data(''' dataname ''', n,noise); % generate ' dataname ' data']};
    
            
            
    delete(hands_l_d.figure1);
    drawnow;
else
    if get(hands_l_d.file,'value')
        % here if need to load data from file
        handles.isxls=false;
        
        delete(hands_l_d.figure1);
        drawnow;

        S = uiimport;
        if length(S)==0
            handles.X=[];
            handles.labels=[];
            handles.islb=false;
            handles.isl=false;
            set(handles.dld,'visible','off');
            handles.m1={}; 
            handles.m2={}; 
            handles.m3={};
            handles.mstf={};
        else
            vn=fieldnames(S);
            if length(vn)==1
                X = getfield(S, vn{1});
                hs1=load_data_1_var;
                hld1=guidata(hs1);
                if get(hld1.i,'value')
                    cn=str2num(get(hld1.cn,'string'));
                    lX=length(X(1,:));
                    ncn1=1:lX;
                    ncni=find(ncn1~=cn);
                    ncn=ncn1(ncni);
                    handles.X=X(:,ncn);
                    handles.isl=true;
                    set(handles.dld,'visible','on');
                    labels=X(:,cn);
                    handles.labels=labels;
                    handles.islb=true;
                else
                    % one variable without labels
                    handles.X=X;
                    handles.isl=true;
                    set(handles.dld,'visible','on');
                    handles.labels=[];
                    handles.islb=false;
                    
                    
                    
                end
                delete(hld1.figure1);
            else
                 hss=load_data_vars(vn);
                 hlds=guidata(hss);
                 Xn=get(hlds.d,'string');
                 X = getfield(S, Xn);
                 if get(hlds.ni,'value') % not include
                     handles.labels=[];
                     handles.islb=false;
                     handles.X=X;
                     handles.isl=true;
                     set(handles.dld,'visible','on');
                 end
                 if get(hlds.iv,'value') % include from variable
                     ivn=get(hlds.ivn,'String');
                     labels = getfield(S, ivn);
                     handles.labels=labels;
                     handles.islb=true;
                     handles.X=X;
                     handles.isl=true;
                     set(handles.dld,'visible','on');
                 end
                 if get(hlds.ic,'value') % include from column
                     cn=str2num(get(hlds.icn,'String'));
                     lX=length(X(1,:));
                     ncn1=1:lX;
                     ncni=find(ncn1~=cn);
                     ncn=ncn1(ncni);
                     handles.X=X(:,ncn);
                     handles.isl=true;
                     set(handles.dld,'visible','on');
                     labels=X(:,cn);
                     handles.labels=labels;
                     handles.islb=true;
                 end
                 delete(hlds.figure1);
            end
            % history:
            % clear:
            handles.m1={}; 
            handles.m2={}; 
            handles.m3={};
            handles.mstf={};
            
            if handles.islb
                handles.m1={'load(''X.mat''); % load data';
                            'load(''labels.mat''); % load labels';};
                         
            else
                handles.m1={'load(''X.mat''); % load data'};
            end
        end
    else
        % here if load from xls file
        handles.isxls=true;
        delete(hands_l_d.figure1);
        drawnow;
        [FileName,PathName] = uigetfile({'*.xls';'*.xlsx'},'Select the xls-file');
        if (FileName==0)
            % do nothing if click cancel
        else
            X = xlsread([PathName FileName]); % load fom xls numbers only
            
            handles.m1={}; 
            handles.m2={}; 
            handles.m3={};
            handles.mstf={};
            
            hstmp=load_xls(length(X(1,:)));
            drawnow;
            hands_l_x=guidata(hstmp);

            if get(hands_l_x.wl,'value')
                % if not use labels
                handles.labels=[];
                handles.islb=false;
                handles.X=X;
                handles.isl=true;
                
                % history:
                handles.m1={['PathName = ''' PathName '''; % path to xls-file'];
                            ['FileName = ''' FileName '''; % file name of xls-file'];
                            'X = xlsread([PathName FileName]); % load xls file'};
                
                set(handles.dld,'visible','on');
            else
                % if use labels
                cn=str2num(get(hands_l_x.col,'String'));
                lX=length(X(1,:));
                ncn1=1:lX;
                ncni=find(ncn1~=cn);
                ncn=ncn1(ncni);
                handles.X=X(:,ncn);
                handles.isl=true;
                set(handles.dld,'visible','on');
                labels=X(:,cn);
                handles.labels=labels;
                handles.islb=true;
                
                % history:
                handles.m1={['PathName = ''' PathName '''; % path to xls-file'];
                            ['FileName = ''' FileName '''; % file name of xls-file'];
                            'X = xlsread([PathName FileName]); % load xls file';
                            ['cn = ' num2str(cn) '; % column number where labels are placed'];
                            'lX=length(X(1,:)); % total number of column';
                            'ncn1=1:lX;';
                            'ncni=find(ncn1~=cn); % indexes of data columns';
                            'ncn=ncn1(ncni); % data columns';
                            'labels=X(:,cn); % get labels';
                            'X=X(:,ncn); % get data'};
            end

            delete(hands_l_x.figure1);
            drawnow;
        end
    end
end

set(handles.edr,'String','');

set(handles.cd,'visible','off');
handles.mcd=false;


guidata(handles.figure1, handles);

set(handles.ld,'Enable','on');
drawnow;


% --- Executes on button press in sp.
function sp_Callback(hObject, eventdata, handles)
if handles.isl
    ld=length(handles.X(1,:));
    if ld==2
        hf=figure;
        set(hf,'name','Original dataset','NumberTitle','off');
        if handles.islb
            scatter(handles.X(:,1),handles.X(:,2),5,handles.labels);
        else
            scatter(handles.X(:,1),handles.X(:,2),5);
        end
        title('Original dataset');
    else

        if handles.islb
            scattern('Original dataset',handles.X,handles.labels);
        else
            scattern('Original dataset',handles.X);
        end
        
    end

else
    % 'not loaded'
    not_loaded;
end

% --- Executes on button press in p.
function p_Callback(hObject, eventdata, handles)
if handles.isl
    ld=length(handles.X(1,:));
    if ld==2
        hf=figure;
        set(hf,'name','Original dataset','NumberTitle','off');
%         if handles.islb
%             plot(handles.X(:,1),handles.X(:,2),5,handles.labels);
%         else
            plot(handles.X(:,1),handles.X(:,2),'.r');
%         end
        title('Original dataset');
    else

%         if handles.islb
%             scattern('Original dataset',handles.X,handles.labels);
%         else
            plotn('Original dataset',handles.X);
%         end
        
    end

else
    % 'not loaded'
    not_loaded;
end


% --- Executes on button press in ed.
function ed_Callback(hObject, eventdata, handles)
handles.a2=false;
if handles.isl
    try
        s=choose_method;
    catch
        return
    end
    hs1=guidata(s);
    set(hs1.ok,'Enable','off');
    drawnow;
    method=get(hs1.str,'string');
    no_dims = intrinsic_dim(handles.X, method);
    
    % estimate dimention lead to clear calculated data:
    handles.mcd=false; % if mapping was calculated
    handles.mX=[]; % mapped X
    set(handles.cd,'Visible','off');
    
    % history:
    handles.m2={};
    handles.m3={};
    handles.mstf={};
    
    % get detailed method name from listbox:
    lstbs=get(hs1.listbox1,'string');
    lstv=get(hs1.listbox1,'value');
    mthds=lstbs{lstv};
    
    handles.m2={['method_ed = ''' method '''; % (' mthds ') method of estimation of dimensionality'];
                ['no_dims = intrinsic_dim(X, method_ed); % estimate intrinsic dimensionality']};
    
    
    delete(hs1.figure1);
    drawnow;
    set(handles.edr,'string',num2str(no_dims));
    set(handles.edr,'UserData',no_dims); % memorize pricise value in userdata
    
    
   guidata(handles.figure1, handles);
else
    % 'not loaded'
    not_loaded;
end



function edr_Callback(hObject, eventdata, handles)
% hObject    handle to edr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edr as text
%        str2double(get(hObject,'String')) returns contents of edr as a double


% --- Executes during object creation, after setting all properties.
function edr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in cm.
function cm_Callback(hObject, eventdata, handles)
if handles.isl
    handles.m3={}; % eraise priviose history
    handles.mstf={};
    
    no_dims=get(handles.edr,'UserData');
    no_dims_old=no_dims;
    
    try
        s=mapping_parameters(no_dims,handles.islb);
    catch
        return
    end
    hs1=guidata(s);
    
    %mappedA = compute_mapping(A, type, no_dims, parameters, eig_impl)
    
    no_dims=str2num(get(hs1.nd,'string'));
    %if (~handles.a2)||(round(no_dims_old)~=no_dims) % if was not added or if changed
        handles.a2=true;
        if ~isempty(no_dims_old)
            if round(no_dims_old)~=no_dims
                % if demetion was changed
                handles.m21={' ';
                    ['no_dims = ' num2str(no_dims) '; % supposed number of dimensions'];
                    ' '};
            else
                handles.m21={' ';
                    ['no_dims = round(no_dims); % round number of dimensions to have integer number'];
                    ' '};
            end
        else
            handles.m21={' ';
                ['no_dims = ' num2str(no_dims) '; % supposed number of dimensions'];
                ' '};
        end
        
        if isempty(handles.m2)
            handles.m21={' ';
                    ['no_dims = ' num2str(no_dims) '; % supposed number of dimensions'];
                    ' '};
        end
        %handles.m2=vertcat(handles.m2,m2t);
    %end
    
    mappedA=[];
  try
    noparam=false; % if no parameters
    mthd=''; % method when no parameters
    switch get(hs1.listbox1,'Value')
        case 1 % PCA
            % no parameters;
            mappedA = compute_mapping(handles.X, 'PCA', no_dims);
            noparam=true;
            mthd='PCA';
        case 2 % LDA
            % no parameters;
            
            % correct lables only to column-vector:
             lb=handles.labels;
             slb=size(lb);
             if min(slb)>1
                 warning('slb must be vector');
             end
             if slb(1)<slb(2)
                 lb=lb';
             end
                 
                 
             if handles.islb
                 mappedA = compute_mapping([lb handles.X], 'LDA', no_dims); % set labels through first column
                 if slb(1)<slb(2)
                    handles.m3={['labels=labels''; % labels must be a vector-column '];
                                ['mappedX = compute_mapping([labels X], ''LDA'', no_dims); % set labels through first column']};
                 else
                    handles.m3={['mappedX = compute_mapping([labels X], ''LDA'', no_dims); % set labels through first column']};
                 end
                     
             else
                 % imposible because data without labels
                 warning('it is imposible to be here');
             end
        case 3 % MDS
            % no parameters;
            mappedA = compute_mapping(handles.X, 'MDS', no_dims);
            noparam=true;
            mthd='MDS';
        case 4 % SimplePCA
            % no parameters;
            mappedA = compute_mapping(handles.X, 'SimplePCA', no_dims);
            noparam=true;
            mthd='SimplePCA';
        case 5 % ProbPCA
         
            mi=str2num(get(hs1.mi,'string'));

            mappedA = compute_mapping(handles.X, 'ProbPCA', no_dims, mi);
            
            handles.m3={['mi = ' num2str(mi) '; % max iterations'];
                ['mappedX = compute_mapping(X, ''ProbPCA'', no_dims, mi);']};
        case 6 % FactorAnalysis
            % no parameters;
            mappedA = compute_mapping(handles.X, 'FactorAnalysis', no_dims);
            noparam=true;
            mthd='FactorAnalysis';
        case 7 % Isomap
            if get(hs1.ka,'value')
                k='adaptive';
                handles.m3={['k = ''adaptive''; % adaptive number of  nearest neighbors in a neighborhood graph']};
            else
                k=str2num(get(hs1.k,'string'));
                
                handles.m3={['k = ' num2str(k) '; % number of  nearest neighbors in a neighborhood graph']};
            end
            
            mappedA = compute_mapping(handles.X, 'Isomap', no_dims, k);
            
            m3t={['mappedX = compute_mapping(X, ''Isomap'', no_dims, k);']};
            handles.m3=vertcat(handles.m3,m3t);
        case 8 % LandmarkIsomap
            if get(hs1.ka,'value')
                k='adaptive';
                handles.m3={['k = ''adaptive''; % adaptive number of  nearest neighbors in a neighborhood graph']};
            else
                k=str2num(get(hs1.k,'string'));
                handles.m3={['k = ' num2str(k) '; % number of  nearest neighbors in a neighborhood graph']};
            end
            
            prc=str2num(get(hs1.prc,'string'));
            
            m3t={['prc = ' num2str(prc) '; % percentage']};
            handles.m3=vertcat(handles.m3,m3t);
            
            mappedA = compute_mapping(handles.X, 'LandmarkIsomap', no_dims, k, prc);
            
            m3t={['mappedX = compute_mapping(X, ''LandmarkIsomap'', no_dims, k, prc);']};
            handles.m3=vertcat(handles.m3,m3t);
        case 9 % LLE
            if get(hs1.ka,'value')
                k='adaptive';
                handles.m3={['k = ''adaptive''; % adaptive number of  nearest neighbors in a neighborhood graph']};
            else
                k=str2num(get(hs1.k,'string'));
                handles.m3={['k = ' num2str(k) '; % number of  nearest neighbors in a neighborhood graph']};
            end
            
            if get(hs1.eim,'value')
                eim='Matlab';
                m3t={'eim = ''Matlab''; % eigenanalysis implementation'};
            else
                eim='JDQR';
                m3t={'eim = ''JDQR''; % eigenanalysis implementation'};
            end
            handles.m3=vertcat(handles.m3,m3t);
            
            [mappedA, mapping] = compute_mapping(handles.X, 'LLE', no_dims, k, eim);
            m3t={['[mappedX, mapping] = compute_mapping(X, ''LLE'', no_dims, k, eim);']};
            handles.m3=vertcat(handles.m3,m3t);
            
        case 10 % Laplacian
            if get(hs1.ka,'value')
                k='adaptive';
                handles.m3={['k = ''adaptive''; % adaptive number of  nearest neighbors in a neighborhood graph']};
            else
                k=str2num(get(hs1.k,'string'));
                handles.m3={['k = ' num2str(k) '; % number of  nearest neighbors in a neighborhood graph']};
            end
            
            sig=str2num(get(hs1.sig,'string'));
            m3t={['sig = ' num2str(sig) '; % variance of a Gaussian kernel']}; 
            handles.m3=vertcat(handles.m3,m3t);
            
            if get(hs1.eim,'value')
                eim='Matlab';
                m3t={'eim = ''Matlab''; % eigenanalysis implementation'};
            else
                eim='JDQR';
                m3t={'eim = ''JDQR''; % eigenanalysis implementation'};
            end
            handles.m3=vertcat(handles.m3,m3t);
            
            mappedA = compute_mapping(handles.X, 'Laplacian', no_dims, k, sig, eim);
            m3t={['mappedX = compute_mapping(X, ''Laplacian'', no_dims, k, sig, eim);']};
            handles.m3=vertcat(handles.m3,m3t);
        case 11 % HessianLLE
            if get(hs1.ka,'value')
                k='adaptive';
                handles.m3={['k = ''adaptive''; % adaptive number of  nearest neighbors in a neighborhood graph']};
            else
                k=str2num(get(hs1.k,'string'));
                handles.m3={['k = ' num2str(k) '; % number of  nearest neighbors in a neighborhood graph']};
            end
            
            if get(hs1.eim,'value')
                eim='Matlab';
                m3t={'eim = ''Matlab''; % eigenanalysis implementation'};
            else
                eim='JDQR';
                m3t={'eim = ''JDQR''; % eigenanalysis implementation'};
            end
            handles.m3=vertcat(handles.m3,m3t);
            
            mappedA = compute_mapping(handles.X, 'HessianLLE', no_dims, k, eim);
            m3t={['mappedX = compute_mapping(X, ''HessianLLE'', no_dims, k, eim);']};
            handles.m3=vertcat(handles.m3,m3t);
        case 12 % LTSA
            if get(hs1.ka,'value')
                k='adaptive';
                handles.m3={['k = ''adaptive''; % adaptive number of  nearest neighbors in a neighborhood graph']};
            else
                k=str2num(get(hs1.k,'string'));
                handles.m3={['k = ' num2str(k) '; % number of  nearest neighbors in a neighborhood graph']};
            end
            
            if get(hs1.eim,'value')
                eim='Matlab';
                m3t={'eim = ''Matlab''; % eigenanalysis implementation'};
            else
                eim='JDQR';
                m3t={'eim = ''JDQR''; % eigenanalysis implementation'};
            end
            handles.m3=vertcat(handles.m3,m3t);
            
            mappedA = compute_mapping(handles.X, 'LTSA', no_dims, k, eim);
            m3t={['mappedX = compute_mapping(X, ''LTSA'', no_dims, k, eim);']};
            handles.m3=vertcat(handles.m3,m3t);
        case 13 % MVU
            if get(hs1.ka,'value')
                k='adaptive';
                handles.m3={['k = ''adaptive''; % adaptive number of  nearest neighbors in a neighborhood graph']};
            else
                k=str2num(get(hs1.k,'string'));
                handles.m3={['k = ' num2str(k) '; % number of  nearest neighbors in a neighborhood graph']};
            end
            
            if get(hs1.eim,'value')
                eim='Matlab';
                m3t={'eim = ''Matlab''; % eigenanalysis implementation'};
            else
                eim='JDQR';
                m3t={'eim = ''JDQR''; % eigenanalysis implementation'};
            end
            handles.m3=vertcat(handles.m3,m3t);
            
            mappedA = compute_mapping(handles.X, 'MVU', no_dims, k, eim);
            m3t={['mappedX = compute_mapping(X, ''MVU'', no_dims, k, eim);']};
            handles.m3=vertcat(handles.m3,m3t);
        case 14 % CCA
            if get(hs1.ka,'value')
                k='adaptive';
                handles.m3={['k = ''adaptive''; % adaptive number of  nearest neighbors in a neighborhood graph']};
            else
                k=str2num(get(hs1.k,'string'));
                handles.m3={['k = ' num2str(k) '; % number of  nearest neighbors in a neighborhood graph']};
            end
            
            if get(hs1.eim,'value')
                eim='Matlab';
                m3t={'eim = ''Matlab''; % eigenanalysis implementation'};
            else
                eim='JDQR';
                m3t={'eim = ''JDQR''; % eigenanalysis implementation'};
            end
            handles.m3=vertcat(handles.m3,m3t);
            
            mappedA = compute_mapping(handles.X, 'CCA', no_dims, k, eim);
            m3t={['mappedX = compute_mapping(X, ''CCA'', no_dims, k, eim);']};
            handles.m3=vertcat(handles.m3,m3t);
        case 15 % LandmarkMVU
            if get(hs1.ka,'value')
                k='adaptive';
                handles.m3={['k = ''adaptive''; % adaptive number of  nearest neighbors in a neighborhood graph']};
            else
                k=str2num(get(hs1.k,'string'));
                handles.m3={['k = ' num2str(k) '; % number of  nearest neighbors in a neighborhood graph']};
            end
            
            mappedA = compute_mapping(handles.X, 'LandmarkMVU', no_dims, k);
            m3t={['mappedX = compute_mapping(X, ''LandmarkMVU'', no_dims, k);']};
            handles.m3=vertcat(handles.m3,m3t);
        case 16 % FastMVU
            if get(hs1.ka,'value')
                k='adaptive';
                handles.m3={['k = ''adaptive''; % adaptive number of  nearest neighbors in a neighborhood graph']};
            else
                k=str2num(get(hs1.k,'string'));
                handles.m3={['k = ' num2str(k) '; % number of  nearest neighbors in a neighborhood graph']};
            end
            
            ft=get(hs1.fty,'value');
            if ft
                m3t={['ft = true; % finetune']};
            else
                m3t={['ft = false; % finetune']};
            end
            handles.m3=vertcat(handles.m3,m3t);
            
            if get(hs1.eim,'value')
                eim='Matlab';
                m3t={'eim = ''Matlab''; % eigenanalysis implementation'};
            else
                eim='JDQR';
                m3t={'eim = ''JDQR''; % eigenanalysis implementation'};
            end
            handles.m3=vertcat(handles.m3,m3t);
            
            mappedA = compute_mapping(handles.X, 'FastMVU', no_dims, k, ft, eim);
            m3t={['mappedX = compute_mapping(X, ''FastMVU'', no_dims, k, ft, eim);']};
            handles.m3=vertcat(handles.m3,m3t);
        case 17 % DiffusionMaps
            
            t=str2num(get(hs1.t,'string'));
            handles.m3={['t = ' num2str(t) ';']};
            
            sig=str2num(get(hs1.sig,'string'));
            m3t={['sig = ' num2str(sig) '; % variance of a Gaussian kernel']}; 
            handles.m3=vertcat(handles.m3,m3t);
            
            mappedA = compute_mapping(handles.X, 'DiffusionMaps', no_dims, t, sig);
            m3t={['mappedX = compute_mapping(X, ''DiffusionMaps'', no_dims, t, sig);']};
            handles.m3=vertcat(handles.m3,m3t);
        case 18 % KernelPCA
            kernel='gauss';
            if get(hs1.kl,'value')
                kernel='linear';
                mappedA = compute_mapping(handles.X, 'KernelPCA', no_dims,kernel);
                handles.m3={['kernel = ''linear'';'];
                            ['mappedX = compute_mapping(X, ''KernelPCA'', no_dims, kernel);']};
            end
            if get(hs1.kg,'value')
                s=str2num(get(hs1.kgs,'string'));
                handles.m3={['s = ' num2str(s) '; % variance']};
                kernel='gauss';
                mappedA = compute_mapping(handles.X, 'KernelPCA', no_dims,kernel,s);
                m3t={['kernel = ''gauss'';'];
                     ['mappedX = compute_mapping(X, ''KernelPCA'', no_dims, kernel, s);']};
                handles.m3=vertcat(handles.m3,m3t);
            end
            if get(hs1.kp,'value')
                R=str2num(get(hs1.kpR,'string'));
                handles.m3={['R = ' num2str(R) '; % additional value']};
                d=str2num(get(hs1.kpd,'string'));
                m3t={['d = ' num2str(d) '; % power number']};
                handles.m3=vertcat(handles.m3,m3t);
                kernel='poly';
                mappedA = compute_mapping(handles.X, 'KernelPCA', no_dims,kernel,R,d);
                m3t={['kernel = ''poly'';'];
                     ['mappedX = compute_mapping(X, ''KernelPCA'', no_dims, kernel, R, d);']};
                handles.m3=vertcat(handles.m3,m3t);
            end
           
        case 19 % GDA
                        
            % correct lables only to column-vector:
             lb=handles.labels;
             slb=size(lb);
             if min(slb)>1
                 warning('slb must be vector');
             end
             if slb(1)<slb(2)
                 lb=lb';
             end
             
             if slb(1)<slb(2)
                 m3t={['labels=labels''; % labels must be a vector-column ']};
             else
                 m3t={};
             end
             
             if handles.islb
                kernel='gauss';
                if get(hs1.kl,'value')
                    kernel='linear';
                    mappedA = compute_mapping(handles.X, 'KernelPCA', no_dims,kernel);
                    handles.m3={['kernel = ''linear'';'];
                                ['mappedX = compute_mapping(X, ''KernelPCA'', no_dims, kernel);']};
                    handles.m3=vertcat(m3t,handles.m3);
                end
                if get(hs1.kg,'value')
                    s=str2num(get(hs1.kgs,'string'));
                    handles.m3={['s = ' num2str(s) '; % variance']};
                    handles.m3=vertcat(m3t,handles.m3);
                    kernel='gauss';
                    mappedA = compute_mapping(handles.X, 'KernelPCA', no_dims,kernel,s);
                    m3t={['kernel = ''gauss'';'];
                         ['mappedX = compute_mapping(X, ''KernelPCA'', no_dims, kernel, s);']};
                    handles.m3=vertcat(handles.m3,m3t);
                end
                if get(hs1.kp,'value')
                    R=str2num(get(hs1.kpR,'string'));
                    handles.m3={['R = ' num2str(R) '; % additional value']};
                    handles.m3=vertcat(m3t,handles.m3);
                    d=str2num(get(hs1.kpd,'string'));
                    m3t={['d = ' num2str(d) '; % power number']};
                    handles.m3=vertcat(handles.m3,m3t);
                    kernel='poly';
                    mappedA = compute_mapping(handles.X, 'KernelPCA', no_dims,kernel,R,d);
                    m3t={['kernel = ''poly'';'];
                         ['mappedX = compute_mapping(X, ''KernelPCA'', no_dims, kernel, R, d);']};
                    handles.m3=vertcat(handles.m3,m3t);
                end
                 
             else
                 % imposible because data without labels
                 warning('it is imposible to be here');
             end
        case 20 % SNE
         
            prp=str2num(get(hs1.prp,'string'));
                        
            handles.m3={['prp = ' num2str(prp) '; % perplexity']};
                        
            mappedA = compute_mapping(handles.X, 'SNE', no_dims, prp);
            m3t={['mappedX = compute_mapping(X, ''SNE'', no_dims, prp);']};
            handles.m3=vertcat(handles.m3,m3t);
        case 21 % SymSNE
            
            prp=str2num(get(hs1.prp,'string'));
            handles.m3={['prp = ' num2str(prp) '; % perplexity']};
            
            mappedA = compute_mapping(handles.X, 'SymSNE', no_dims, prp);
            m3t={['mappedX = compute_mapping(X, ''SymSNE'', no_dims, prp);']};
            handles.m3=vertcat(handles.m3,m3t);
        case 22 % t-SNE
            
            prp=str2num(get(hs1.prp,'string'));
            handles.m3={['prp = ' num2str(prp) '; % perplexity']};
            
            mappedA = compute_mapping(handles.X, 't-SNE', no_dims, prp);
            m3t={['mappedX = compute_mapping(X, ''t-SNE'', no_dims, prp);']};
            handles.m3=vertcat(handles.m3,m3t);
        case 23 % LPP
            if get(hs1.ka,'value')
                k='adaptive';
                handles.m3={['k = ''adaptive''; % adaptive number of  nearest neighbors in a neighborhood graph']};
            else
                k=str2num(get(hs1.k,'string'));
                handles.m3={['k = ' num2str(k) '; % number of  nearest neighbors in a neighborhood graph']};
            end
            
            sig=str2num(get(hs1.sig,'string'));
            m3t={['sig = ' num2str(sig) '; % variance of a Gaussian kernel']}; 
            handles.m3=vertcat(handles.m3,m3t);
            
            if get(hs1.eim,'value')
                eim='Matlab';
                m3t={'eim = ''Matlab''; % eigenanalysis implementation'};
            else
                eim='JDQR';
                m3t={'eim = ''JDQR''; % eigenanalysis implementation'};
            end
            handles.m3=vertcat(handles.m3,m3t);
            
            mappedA = compute_mapping(handles.X, 'LPP', no_dims, k, sig, eim);
            m3t={['mappedX = compute_mapping(X, ''LPP'', no_dims, k, sig, eim);']};
            handles.m3=vertcat(handles.m3,m3t);
        case 24 % NPE
            if get(hs1.ka,'value')
                k='adaptive';
                handles.m3={['k = ''adaptive''; % adaptive number of  nearest neighbors in a neighborhood graph']};
            else
                k=str2num(get(hs1.k,'string'));
                handles.m3={['k = ' num2str(k) '; % number of  nearest neighbors in a neighborhood graph']};
            end
            
                        
            if get(hs1.eim,'value')
                eim='Matlab';
                m3t={'eim = ''Matlab''; % eigenanalysis implementation'};
            else
                eim='JDQR';
                m3t={'eim = ''JDQR''; % eigenanalysis implementation'};
            end
            handles.m3=vertcat(handles.m3,m3t);
            
            mappedA = compute_mapping(handles.X, 'NPE', no_dims, k, eim);
            m3t={['mappedX = compute_mapping(X, ''NPE'', no_dims, k, eim);']};
            handles.m3=vertcat(handles.m3,m3t);
        case 25 % LLTSA
            if get(hs1.ka,'value')
                k='adaptive';
                handles.m3={['k = ''adaptive''; % adaptive number of  nearest neighbors in a neighborhood graph']};
            else
                k=str2num(get(hs1.k,'string'));
                handles.m3={['k = ' num2str(k) '; % number of  nearest neighbors in a neighborhood graph']};
            end
            
                        
            if get(hs1.eim,'value')
                eim='Matlab';
                m3t={'eim = ''Matlab''; % eigenanalysis implementation'};
            else
                eim='JDQR';
                m3t={'eim = ''JDQR''; % eigenanalysis implementation'};
            end
            handles.m3=vertcat(handles.m3,m3t);
            
            mappedA = compute_mapping(handles.X, 'LLTSA', no_dims, k, eim);
            m3t={['mappedX = compute_mapping(X, ''LLTSA'', no_dims, k, eim);']};
            handles.m3=vertcat(handles.m3,m3t);
        case 26 % SPE
            tp='Global';
            if get(hs1.tg,'value')
                tp='Global';
                mappedA = compute_mapping(handles.X, 'SPE', no_dims, tp);
                handles.m3={'tp = ''Global''; % type of stress function that minimized';
                            ['mappedX = compute_mapping(X, ''SPE'', no_dims, tp);']};
            end
            if get(hs1.tl,'value')
                tp='Local';
                k=str2num(get(hs1.k,'string'));
                mappedA = compute_mapping(handles.X, 'SPE', no_dims, tp, k);
                handles.m3={['k = ' num2str(k) '; % number of nearest neighbors in a neighborhood graph'];
                            'tp = ''Local''; % type of stress function that minimized';
                            ['mappedX = compute_mapping(X, ''SPE'', no_dims, tp, k);']};
            end
            
        case 27 % AutoEncoderRBM
            % no parameters;
            mappedA = compute_mapping(handles.X, 'AutoEncoderRBM', no_dims);
            noparam=true;
            mthd='AutoEncoderRBM';
        case 28 % AutoEncoderEA
            % no parameters;
            mappedA = compute_mapping(handles.X, 'AutoEncoderEA', no_dims);
            noparam=true;
            mthd='AutoEncoderEA';
        case 29 % LLC
            if get(hs1.ka,'value')
                k='adaptive';
                handles.m3={['k = ''adaptive''; % adaptive number of  nearest neighbors in a neighborhood graph']};
            else
                k=str2num(get(hs1.k,'string'));
                handles.m3={['k = ' num2str(k) '; % number of  nearest neighbors in a neighborhood graph']};
            end
            
            na=str2num(get(hs1.na,'string'));
            m3t={['na = ' num2str(na) '; % number of factor analyzers']};
            handles.m3=vertcat(handles.m3,m3t);
            
            mi=str2num(get(hs1.mi,'string'));
            m3t={['mi = ' num2str(mi) '; % max iterations']};
            handles.m3=vertcat(handles.m3,m3t);
            
            if get(hs1.eim,'value')
                eim='Matlab';
                m3t={'eim = ''Matlab''; % eigenanalysis implementation'};
            else
                eim='JDQR';
                m3t={'eim = ''JDQR''; % eigenanalysis implementation'};
            end
            
            mappedA = compute_mapping(handles.X, 'LLC', no_dims, k, na, mi, eim);
            m3t={['mappedX = compute_mapping(X, ''LLC'', no_dims, k, na, mi, eim);']};
            handles.m3=vertcat(handles.m3,m3t);
        case 30 % ManifoldChart
            
            na=str2num(get(hs1.na,'string'));
            handles.m3={['na = ' num2str(na) '; % number of factor analyzers']};
            
            
            mi=str2num(get(hs1.mi,'string'));
            m3t={['mi = ' num2str(mi) '; % max iterations']};
            handles.m3=vertcat(handles.m3,m3t);
            
            if get(hs1.eim,'value')
                eim='Matlab';
                m3t={'eim = ''Matlab''; % eigenanalysis implementation'};
            else
                eim='JDQR';
                m3t={'eim = ''JDQR''; % eigenanalysis implementation'};
            end
            
            mappedA = compute_mapping(handles.X, 'ManifoldChart', no_dims, na, mi, eim);
            m3t={['mappedX = compute_mapping(X, ''ManifoldChart'', no_dims, na, mi, eim);']};
            handles.m3=vertcat(handles.m3,m3t);
        case 31 % CFA
                        
            na=str2num(get(hs1.na,'string'));
            handles.m3={['na = ' num2str(na) '; % number of factor analyzers']};
            
            
            mi=str2num(get(hs1.mi,'string'));
            m3t={['mi = ' num2str(mi) '; % max iterations']};
            handles.m3=vertcat(handles.m3,m3t);
            
            mappedA = compute_mapping(handles.X, 'CFA', no_dims, na, mi);
            m3t={['mappedX = compute_mapping(X, ''CFA'', no_dims, na, mi);']};
            handles.m3=vertcat(handles.m3,m3t);
        case 32 % GPLVM
            sig=str2num(get(hs1.sig,'string'));
            m3t={['sig = ' num2str(sig) '; % variance of a Gaussian kernel']}; 
            handles.m3=vertcat(handles.m3,m3t);
            
            mappedA = compute_mapping(handles.X, 'GPLVM', no_dims, sig);
            m3t={['mappedX = compute_mapping(X, ''GPLVM'', no_dims, sig);']};
            handles.m3=vertcat(handles.m3,m3t);
        case 33 % NCA
            mappedA = compute_mapping(handles.X, 'NCA', no_dims);
            m3t={'mappedX = compute_mapping(X, ''NCA'', no_dims);'};
            handles.m3=vertcat(handles.m3,m3t);
        case 34 % MCML
            mappedA = compute_mapping(handles.X, 'MCML', no_dims);
            m3t={'mappedX = compute_mapping(X, ''MCML'', no_dims);'};
            handles.m3=vertcat(handles.m3,m3t);

    end
    if noparam
        handles.m3={['mappedX = compute_mapping(X, ''' mthd ''', no_dims); % compute mapping using ' mthd ' method']};
    end
  catch
      set(handles.cd,'visible','off');
      handles.mcd=false;
      warning('mapping was not calculated');
      handles.m3={};
      handles.mstf={};
      guidata(handles.figure1, handles);
      delete(hs1.figure1);
      drawnow;
      rethrow(lasterror);
      return
  end
    if length(mappedA)~=0
        set(handles.cd,'visible','on');
        handles.mX=mappedA;
        handles.mcd=true;
    else
        set(handles.cd,'visible','off');
        handles.mcd=false;
        warning('mapping was not calculated');
        handles.m3={};
        handles.mstf={};
    end
    guidata(handles.figure1, handles);
    delete(hs1.figure1);
    drawnow;
    
else
    % 'not loaded'
    not_loaded;
end


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
if handles.mcd
    ld=length(handles.mX(1,:));
    if ld==2
        hf=figure;
        set(hf,'name','Result of dimensionality reduction','NumberTitle','off');
        if handles.islb && size(handles.mX, 1) == numel(handles.labels)
            scatter(handles.mX(:,1),handles.mX(:,2),5,handles.labels);
        else
            scatter(handles.mX(:,1),handles.mX(:,2),5);
        end
        title('Result of dimensionality reduction');
        if size(handles.mX, 1) ~= numel(handles.labels)
            warning('The GUI cannot yet deal properly with disconnected parts in the neighborhood graph.');
        end
    else
        if ld==1
            hf=figure;
            set(hf,'name','Result of dimensionality reduction','NumberTitle','off');
            if handles.islb && size(handles.mX, 1) == numel(handles.labels)
                scatter(handles.mX(:,1),zeros(length(handles.mX(:,1)),1),5,handles.labels);
            else
                scatter(handles.mX(:,1),zeros(length(handles.mX(:,1)),1),5);
            end
            title('Result of dimensionality reduction');
            if size(handles.mX, 1) ~= numel(handles.labels)
                warning('The GUI cannot yet deal properly with disconnected parts in the neighborhood graph.');
            end
        else
            if handles.islb
                scattern('Result of dimensionality reduction',handles.mX,handles.labels);
            else
                scattern('Result of dimensionality reduction',handles.mX);
            end
        end
        
    end

else
    % 'not calulated'
    not_calculated;
end


% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
if handles.mcd
    ld=length(handles.mX(1,:));
    if ld==2
        hf=figure;
        set(hf,'name','Result of dimensionality reduction','NumberTitle','off');
%         if handles.islb
%             scatter(handles.mX(:,1),handles.mX(:,2),5,handles.labels);
%         else
            plot(handles.mX(:,1),handles.mX(:,2),'.r');
%         end
        title('Result of dimensionality reduction');
    else
        if ld==1
            hf=figure;
            set(hf,'name','Result of dimensionality reduction','NumberTitle','off');
            %if handles.islb
                %scatter(handles.mX(:,1),zeros(length(handles.mX(:,1)),1),5,handles.labels);
            %else
                plot(handles.mX(:,1),zeros(length(handles.mX(:,1)),1),'.k');
            %end
            title('Result of dimensionality reduction');
        else
            %if handles.islb
               % scattern('Result of dimensionality reduction',handles.mX,handles.labels);
            %else
                plotn('Result of dimensionality reduction',handles.mX);
            %end
        end
        
    end

else
    % 'not calulated'
    not_calculated;
end


% --- Executes on button press in stmf.
function stmf_Callback(hObject, eventdata, handles)
[file,path] = uiputfile('*.mat','Save Mapped Data As');
if length(file)==1
    if file==0
        return
    end
end
mX=handles.mX;
save([path file], 'mX');


% --- Executes on button press in sttf.
function sttf_Callback(hObject, eventdata, handles)
[file,path] = uiputfile('*.txt','Save Mapped Data As');
if length(file)==1
    if file==0
        return
    end
end
mX=handles.mX;
save([path file], 'mX','-ascii', '-tabs');




% --- Executes on button press in stf.
function stf_Callback(hObject, eventdata, handles)
set(handles.stf,'Enable','off');
drawnow;
if handles.mcd
    if get(handles.stmfr,'value')
        [file,path] = uiputfile('*.mat','Save Mapped Data As');
        if length(file)==1
            if file==0
                set(handles.stf,'Enable','on');
                drawnow;
                return
            end
        end
        mX=handles.mX;
        save([path file], 'mX');
        
        handles.mstf={['save(''' path file ''', ''mappedX''); % save result to mat-file']};
    end

    if get(handles.sttfr,'value')
        [file,path] = uiputfile('*.txt','Save Mapped Data As');
        if length(file)==1
            if file==0
                set(handles.stf,'Enable','on');
                drawnow;
                return
            end
        end
        mX=handles.mX;
        save([path file], 'mX','-ascii', '-tabs');
        handles.mstf={['save(''' path file ''', ''mappedX'',''-ascii'', ''-tabs''); % save result to txt-file']};
    end

    if get(handles.stxfr,'value')
        [file,path] = uiputfile('*.xls','Save Mapped Data As');
        if length(file)==1
            if file==0
                set(handles.stf,'Enable','on');
                drawnow;
                return
            end
        end
        mX=handles.mX;
        xlswrite([path file], mX);
        
        handles.mstf={['xlswrite(''' path file ''', mappedX); % save result to xls-file']};

    end
    
    guidata(handles.figure1, handles);
else
    % 'not calulated'
    not_calculated;
end

set(handles.stf,'Enable','on');
drawnow;
    
% --- Executes when selected object is changed in uipanel1.
function uipanel1_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel1 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)




% --- Executes on button press in sh.
function sh_Callback(hObject, eventdata, handles)

set(handles.sh,'enable','off');
drawnow;

if (isempty(handles.m1))&&(isempty(handles.m2))&&(isempty(handles.m3))
    
    no_history;
    
else

    [file,path] = uiputfile('*.m','Save History As');
    
    if length(file)==1
        if file==0
            set(handles.sh,'enable','on');
            drawnow;
            return
        end
    end

    fid = fopen([path file],'w');
    
    fprintf(fid,'%s\r\n',['% this file was generated by drtool at ' datestr(now)]);
    
    fprintf(fid,'%s\r\n',' '); % empty line - delimeter
    
    if ~isempty(handles.isxls)
        if ~handles.isxls
            % here if txt/mat file was loaded so it is need save data and labels
            X=handles.X;
            save([path 'X.mat'],'X');
            if handles.islb 
                % save labells if any
                labels=handles.labels;
                save([path 'labels.mat'],'labels');
                fprintf(fid,'%s\r\n','% Note: data and labels were saved in X.mat and labels.mat together with this m-file in the same folder');
            else
                fprintf(fid,'%s\r\n','% Note: data was saved in X.mat together with this m-file in the same folder');
            end
        end
    end
    
    fprintf(fid,'%s\r\n',' ');
    fprintf(fid,'%s\r\n','% 1.');
    fprintf(fid,'%s\r\n','% get data');
    for mc=1:length(handles.m1)
        fprintf(fid,'%s\r\n',handles.m1{mc});
    end
    
    % plot original data:
    mpod={'% plot original data:'};
    ld=length(handles.X(1,:));
    if ld==2
        mpodt={'hf=figure;';
              'set(hf,''name'',''Original dataset'',''NumberTitle'',''off'');'};
        mpod=vertcat(mpod,mpodt);
        if handles.islb
            mpodt={'scatter(X(:,1),X(:,2),5,labels);'};
        else
            mpodt={'plot(X(:,1),X(:,2),''x-'');'};
        end
        mpod=vertcat(mpod,mpodt);
        
        mpodt={'title(''Original dataset'');'};
        mpod=vertcat(mpod,mpodt);
    else

        if handles.islb
            mpodt={'scattern(''Original dataset'',X,labels);'};
        else
            mpodt={'plotn(''Original dataset'',X);'};
        end
        mpod=vertcat(mpod,mpodt);
    end
    
    fprintf(fid,'%s\r\n',' ');
    for mc=1:length(mpod)
        fprintf(fid,'%s\r\n',mpod{mc});
    end
    fprintf(fid,'%s\r\n',' ');


    
    
    if length(handles.m2)~=0
        fprintf(fid,'%s\r\n',' '); % empty line - delimeter
        
        handles.m2=vertcat(handles.m2,handles.m21);
        
        fprintf(fid,'%s\r\n','% 2.');
        fprintf(fid,'%s\r\n','% estimate intrinsic dimensionality');
        for mc=1:length(handles.m2)
            fprintf(fid,'%s\r\n',handles.m2{mc});
        end
    else
        if length(handles.m3)~=0
            % if was not dimetion estimation thet it is need to set it for
            % part 3
            for mc=1:length(handles.m21)
                fprintf(fid,'%s\r\n',handles.m21{mc});
            end
        end
    end
    
    if length(handles.m3)~=0
    
        fprintf(fid,'%s\r\n',' '); % empty line - delimeter


        fprintf(fid,'%s\r\n','% 3.');
        fprintf(fid,'%s\r\n','% compute mapping');
        for mc=1:length(handles.m3)
            fprintf(fid,'%s\r\n',handles.m3{mc});
        end
        
        % plot result
        mpod={'% plot result of dimensionality reduction:'};
        if handles.islb
            mpodt={'scatter12n(''Result of dimensionality reduction'',mappedX,labels);'};
        else
            mpodt={'plot12n(''Result of dimensionality reduction'',mappedX);'};
        end
        mpod=vertcat(mpod,mpodt);
        fprintf(fid,'%s\r\n',' ');
        for mc=1:length(mpod)
            fprintf(fid,'%s\r\n',mpod{mc});
        end
        fprintf(fid,'%s\r\n',' ');
    end
    
    if length(handles.mstf)~=0
    
        fprintf(fid,'%s\r\n',' '); % empty line - delimeter
        fprintf(fid,'%s\r\n',' '); % empty line - delimeter


        
        for mc=1:length(handles.mstf)
            fprintf(fid,'%s\r\n',handles.mstf{mc});
        end
    end

    fclose(fid);
    
end

set(handles.sh,'enable','on');
drawnow;
