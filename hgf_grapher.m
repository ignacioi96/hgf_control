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

% Last Modified by GUIDE v2.5 13-Jun-2018 12:31:47

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


init_values(handles);

handles.func = varargin{1};


if(get(handles.calc_box, 'Value'))     
    plotter(handles); 
end

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

function val = edit_to_double(edit_box)
val = str2double(edit_box.get('String'));

function val = funky(x)
val = sin(x)+cos(5*x);

function val = linear(x)
val = x;

function init_values(handles)
set(handles.time_edit, 'String', num2str(500.));
set(handles.belief_lambda_edit, 'String', num2str(0));
set(handles.belief_alpha_edit, 'String', num2str(1));
set(handles.belief_omega_edit, 'String', num2str(0.5));
set(handles.belief_kappa_edit, 'String', num2str(0.5));
set(handles.belief_theta_edit, 'String', num2str(1));
set(handles.actual_lambda_edit, 'String', num2str(0));
set(handles.actual_alpha_edit, 'String', num2str(1));
set(handles.mu_des_edit, 'String', num2str(0));
set(handles.pi_des_edit, 'String', num2str(0.1));
set(handles.env_effect_edit, 'String', num2str(0.05));
set(handles.x_init_edit, 'String', num2str(5));
set(handles.mu1_init_edit, 'String', num2str(1));
set(handles.mu2_init_edit, 'String', num2str(1));
set(handles.mu_init_gaussian_edit, 'String', num2str(1));
set(handles.pi_init_gaussian_edit, 'String', num2str(.5));
set(handles.env_effect_period_edit, 'String', num2str(.1));


%% Plotting function
function plotter(handles)
time_val = edit_to_double(handles.time_edit);
belief_lambda_val = edit_to_double(handles.belief_lambda_edit);
belief_alpha_val = edit_to_double(handles.belief_alpha_edit);
belief_omega_val = edit_to_double(handles.belief_omega_edit);
belief_kappa_val = edit_to_double(handles.belief_kappa_edit);
actual_lambda_val = edit_to_double(handles.actual_lambda_edit);
actual_alpha_val = edit_to_double(handles.actual_alpha_edit);
belief_theta_val = edit_to_double(handles.belief_theta_edit);
env_effect_val = edit_to_double(handles.env_effect_edit);
mu_des_val = edit_to_double(handles.mu_des_edit);
pi_des_val = edit_to_double(handles.pi_des_edit);
x_init_val = edit_to_double(handles.x_init_edit);
mu1_init_val = edit_to_double(handles.mu1_init_edit);
mu2_init_val = edit_to_double(handles.mu2_init_edit);
mu_init_gaussian = edit_to_double(handles.mu_init_gaussian_edit);
pi_init_gaussian = edit_to_double(handles.pi_init_gaussian_edit);
env_effect_period_val = edit_to_double(handles.env_effect_period_edit);

str = handles.env_effect_functions.get('String');
val = handles.env_effect_functions.get('Value');
cell_str = str(val);

switch cell_str{1}
    case 'sin(t)'
        env_effect_func = @sin;
    case 'exp(t)'
        env_effect_func = @exp;
    case 'log(t)'
        env_effect_func = @log;
    case 'sin(t) + cos(5t)'
        env_effect_func = @funky;
    case 'linear(t)'
        env_effect_func = @linear;
    case 'tanh(t)'
        env_effect_func = @tanh;
end

action_options = handles.action_popup.get('String');
action_type = action_options(handles.action_popup.get('Value'));
model_options = handles.model_popup.get('String');
model_type = model_options(handles.model_popup.get('Value'));

set(handles.S_prediction_sum_text_old, 'String',...
    num2str(edit_to_double(handles.S_prediction_sum_text)));
set(handles.S_control_sum_text_old, 'String',...
    num2str(edit_to_double(handles.S_control_sum_text)));
set(handles.mean_sq_error_sum_text_old, 'String',...
    num2str(edit_to_double(handles.mean_sq_error_sum_text)));

[u, mus, x, actions, env_effects, action_effects,...
    S_prediction, S_control, mean_sq_error] = handles.func(time_val,...
    belief_lambda_val,belief_alpha_val, belief_omega_val,...
    belief_kappa_val,actual_lambda_val, actual_alpha_val,...
    belief_theta_val, env_effect_val,mu_des_val, pi_des_val, x_init_val,...
    mu1_init_val, mu2_init_val, mu_init_gaussian, pi_init_gaussian,...
    env_effect_func, model_type, action_type, env_effect_period_val);

set(handles.S_prediction_sum_text, 'String',...
    num2str(round(sum(S_prediction),2)));
set(handles.S_control_sum_text, 'String',...
    num2str(round(sum(S_control),2)));
set(handles.mean_sq_error_sum_text, 'String',...
    num2str(round(sum(mean_sq_error),2)));


axes(handles.axes1);
cla;
plot(x, 'b');
hold on;
plot(mus(1,:), 'r');
hold on;
plot(u, 'g');
hold on;
xlim([0 time_val+1]);
title('System State');
axis square;
legend('Real Value (x)', 'Believed Value (mu1)', 'Perceived Value (u)');%,...
%     'Value of mean volatility (mu2)');

axes(handles.axes2);
cla;
plot(actions);
hold on;
plot(action_effects, 'r');
hold on;
plot(env_effects, 'g');
xlim([0 time_val+1]);
legend('Agent Actions', 'Action Effects', 'Env. Perturbations');
title('Turn-by-turn changes');
axis square;

%% Callback and Create Functions
% --- Executes on slider movement.
function belief_lambda_edit_Callback(hObject, eventdata, handles)
% hObject    handle to belief_lambda_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
set(handles.belief_lambda_edit, 'String',(get(hObject,'String')));
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
set(handles.belief_alpha_edit,'String', (get(hObject,'String')));
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
set(handles.belief_omega_edit, 'String',(get(hObject,'String')));
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
set(handles.actual_lambda_edit,'String',(get(hObject,'String')));
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
set(handles.actual_alpha_edit,'String',(get(hObject,'String')));
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
set(handles.belief_kappa_edit,'String',(get(hObject,'String')));

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



function mu_des_edit_Callback(hObject, eventdata, handles)
% hObject    handle to mu_des_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mu_des_edit as text
%        str2double(get(hObject,'String')) returns contents of mu_des_edit as a double
set(handles.mu_des_edit,'String',(get(hObject,'String')));

plotter(handles);

% --- Executes during object creation, after setting all properties.
function mu_des_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mu_des_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function pi_des_edit_Callback(hObject, eventdata, handles)
% hObject    handle to pi_des_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pi_des_edit as text
%        str2double(get(hObject,'String')) returns contents of pi_des_edit as a double
set(handles.pi_des_edit, 'String',(get(hObject,'String')));

plotter(handles);

% --- Executes during object creation, after setting all properties.
function pi_des_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pi_des_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function env_effect_edit_Callback(hObject, eventdata, handles)
% hObject    handle to env_effect_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of env_effect_edit as text
%        str2double(get(hObject,'String')) returns contents of env_effect_edit as a double
set(handles.env_effect_edit,'String', (get(hObject,'String')));
plotter(handles);

% --- Executes during object creation, after setting all properties.
function env_effect_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to env_effect_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function belief_theta_edit_Callback(hObject, eventdata, handles)
% hObject    handle to belief_theta_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of belief_theta_edit as text
%        str2double(get(hObject,'String')) returns contents of belief_theta_edit as a double
set(handles.belief_theta_edit,'String', (get(hObject,'String')));

plotter(handles);

% --- Executes during object creation, after setting all properties.
function belief_theta_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to belief_theta_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function x_init_edit_Callback(hObject, eventdata, handles)
% hObject    handle to x_init_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of x_init_edit as text
%        str2double(get(hObject,'String')) returns contents of x_init_edit as a double
set(handles.x_init_edit,'String', (get(hObject,'String')));
plotter(handles);

% --- Executes during object creation, after setting all properties.
function x_init_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to x_init_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in reset_val_button.
function reset_val_button_Callback(hObject, eventdata, handles)
% hObject    handle to reset_val_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
init_values(handles);
plotter(handles);



function mu1_init_edit_Callback(hObject, eventdata, handles)
% hObject    handle to mu1_init_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mu1_init_edit as text
%        str2double(get(hObject,'String')) returns contents of mu1_init_edit as a double
set(handles.mu1_init_edit,'String', (get(hObject,'String')));
plotter(handles);

% --- Executes during object creation, after setting all properties.
function mu1_init_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mu1_init_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function mu2_init_edit_Callback(hObject, eventdata, handles)
% hObject    handle to mu2_init_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mu2_init_edit as text
%        str2double(get(hObject,'String')) returns contents of mu2_init_edit as a double
set(handles.mu2_init_edit,'String', (get(hObject,'String')));
plotter(handles);

% --- Executes during object creation, after setting all properties.
function mu2_init_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mu2_init_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in env_effect_functions.
function env_effect_functions_Callback(hObject, eventdata, handles)
% hObject    handle to env_effect_functions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns env_effect_functions contents as cell array
%        contents{get(hObject,'Value')} returns selected item from env_effect_functions
str = get(hObject, 'String');
val = get(hObject, 'Value');
set(handles.env_effect_functions, 'String', str);
set(handles.env_effect_functions, 'Value', val);
plotter(handles);
        
% --- Executes during object creation, after setting all properties.
function env_effect_functions_CreateFcn(hObject, eventdata, handles)
% hObject    handle to env_effect_functions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function mu_init_gaussian_edit_Callback(hObject, eventdata, handles)
% hObject    handle to mu_init_gaussian_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mu_init_gaussian_edit as text
%        str2double(get(hObject,'String')) returns contents of mu_init_gaussian_edit as a double
set(handles.mu_init_gaussian_edit,'String', (get(hObject,'String')));
if(get(handles.calc_box, 'Value'))
    plotter(handles);
end

% --- Executes during object creation, after setting all properties.
function mu_init_gaussian_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mu_init_gaussian_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function pi_init_gaussian_edit_Callback(hObject, eventdata, handles)
% hObject    handle to pi_init_gaussian_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pi_init_gaussian_edit as text
%        str2double(get(hObject,'String')) returns contents of pi_init_gaussian_edit as a double
set(handles.pi_init_gaussian_edit,'String', (get(hObject,'String')));
if(get(handles.calc_box, 'Value'))
    plotter(handles);
end

% --- Executes during object creation, after setting all properties.
function pi_init_gaussian_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pi_init_gaussian_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in model_popup.
function model_popup_Callback(hObject, eventdata, handles)
% hObject    handle to model_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns model_popup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from model_popup
str = get(hObject, 'String');
val = get(hObject, 'Value');
set(handles.model_popup, 'String', str);
set(handles.model_popup, 'Value', val);
if(get(handles.calc_box, 'Value'))
    plotter(handles);
end


% --- Executes during object creation, after setting all properties.
function model_popup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to model_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in action_popup.
function action_popup_Callback(hObject, eventdata, handles)
% hObject    handle to action_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns action_popup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from action_popup
str = get(hObject, 'String');
val = get(hObject, 'Value');
set(handles.action_popup, 'String', str);
set(handles.action_popup, 'Value', val);
if(get(handles.calc_box, 'Value'))
    plotter(handles);
end

% --- Executes during object creation, after setting all properties.
function action_popup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to action_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in calc_button.
function calc_button_Callback(hObject, eventdata, handles)
% hObject    handle to calc_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
plotter(handles);

% --- Executes on button press in calc_box.
function calc_box_Callback(hObject, eventdata, handles)
% hObject    handle to calc_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of calc_box
if(get(handles.calc_box, 'Value'))
    plotter(handles);
end

function env_effect_period_edit_Callback(hObject, eventdata, handles)
% hObject    handle to env_effect_period_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of env_effect_period_edit as text
%        str2double(get(hObject,'String')) returns contents of env_effect_period_edit as a double
set(handles.env_effect_period_edit,'String', (get(hObject,'String')));
if(get(handles.calc_box, 'Value'))
    plotter(handles);
end

% --- Executes during object creation, after setting all properties.
function env_effect_period_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to env_effect_period_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
