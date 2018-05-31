function varargout = hgf_grapher(varargin)
% HGF_GRAPHER MATLAB code for hgf_grapher.fig
%      HGF_GRAPHER, by itself, creates a new HGF_GRAPHER or raises the existing
%      singleton*.
%
%      H = HGF_GRAPHER returns the handle to a new HGF_GRAPHER or the handle to
%      the existing singleton*.
%
%      HGF_GRAPHER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in HGF_GRAPHER.M with the given input arguments.
%
%      HGF_GRAPHER('Property','Value',...) creates a new HGF_GRAPHER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before hgf_grapher_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to hgf_grapher_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help hgf_grapher

% Last Modified by GUIDE v2.5 31-May-2018 20:53:57

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @hgf_grapher_OpeningFcn, ...
                   'gui_OutputFcn',  @hgf_grapher_OutputFcn, ...
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


% --- Executes just before hgf_grapher is made visible.
function hgf_grapher_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to hgf_grapher (see VARARGIN)

handles.time_interval = 250.;
handles.belief_lambda_val = 1.;
handles.belief_alpha_val = 1.;
handles.belief_omega_val = 1.;
handles.belief_kappa_val = 1.;
handles.actual_lambda_val = 1.;
handles.actual_alpha_val = 1.;

handles.func = varargin{1};

plotter(handles);

% Choose default command line output for hgf_grapher
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes hgf_grapher wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = hgf_grapher_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function belief_lambda_edit_Callback(hObject, eventdata, handles)
% hObject    handle to belief_lambda_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.belief_lambda_val = str2double(get(hObject,'String'));
plotter(handles);

% --- Executes during object creation, after setting all properties.
function belief_lambda_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to belief_lambda_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function belief_alpha_edit_Callback(hObject, eventdata, handles)
% hObject    handle to belief_alpha_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: str2double(get(hObject,'String')) returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.belief_alpha_val = str2double(get(hObject,'String'));
plotter(handles);

% --- Executes during object creation, after setting all properties.
function belief_alpha_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to belief_alpha_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function belief_omega_edit_Callback(hObject, eventdata, handles)
% hObject    handle to belief_omega_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.belief_omega_val = str2double(get(hObject,'String'));
plotter(handles);

% --- Executes during object creation, after setting all properties.
function belief_omega_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to belief_omega_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function actual_lambda_edit_Callback(hObject, eventdata, handles)
% hObject    handle to belief_lambda_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.actual_lambda_val = str2double(get(hObject,'String'));
plotter(handles);

% --- Executes during object creation, after setting all properties.
function actual_lambda_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to belief_lambda_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function actual_alpha_edit_Callback(hObject, eventdata, handles)
% hObject    handle to actual_alpha_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.actual_alpha_val = str2double(get(hObject,'String'));
plotter(handles);

% --- Executes during object creation, after setting all properties.
function actual_alpha_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to actual_alpha_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes1


% --- Executes on mouse press over axes background.
function axes1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function axes2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes1

% --- Executes during object creation, after setting all properties.
function var_panel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to var_panel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


function plotter(handles)
[u, mus, x, actions] = handles.func(handles.time_interval,...
    handles.belief_lambda_val,handles.belief_alpha_val,...
    handles.belief_omega_val, handles.belief_kappa_val,...
    handles.actual_lambda_val, handles.actual_alpha_val);

axes(handles.axes1);
cla;
plot(x, 'b');
hold on;
plot(mus(1,:), 'r');
hold on;
plot(u, 'g');
hold on;
title('Values of X');
axis square;
legend('Real Value (x)', 'Believed Value (mu)', 'Perceived Value (u)');

axes(handles.axes2);
plot(actions);
title('Actions');
axis square;



function time_edit_Callback(hObject, eventdata, handles)
% hObject    handle to time_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of time_edit as text
%        str2double(get(hObject,'String')) returns contents of time_edit as a double
handles.time_interval = str2double(get(hObject,'String'));
plotter(handles);

% --- Executes during object creation, after setting all properties.
function time_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to time_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function belief_kappa_edit_Callback(hObject, eventdata, handles)
% hObject    handle to belief_kappa_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of belief_kappa_edit as text
%        str2double(get(hObject,'String')) returns contents of belief_kappa_edit as a double
handles.belief_kappa_val = str2double(get(hObject,'String'));
plotter(handles);

% --- Executes during object creation, after setting all properties.
function belief_kappa_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to belief_kappa_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
