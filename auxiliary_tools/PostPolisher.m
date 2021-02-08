function varargout = PostPolisher(varargin)
% POSTPOLISHER MATLAB code for PostPolisher.fig
%      POSTPOLISHER, by itself, creates a new POSTPOLISHER or raises the existing
%      singleton*.
%
%      H = POSTPOLISHER returns the handle to a new POSTPOLISHER or the handle to
%      the existing singleton*.
%
%      POSTPOLISHER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in POSTPOLISHER.M with the given input arguments.
%
%      POSTPOLISHER('Property','Value',...) creates a new POSTPOLISHER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before PostPolisher_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to PostPolisher_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help PostPolisher

% Last Modified by GUIDE v2.5 08-Feb-2021 17:01:45

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @PostPolisher_OpeningFcn, ...
                   'gui_OutputFcn',  @PostPolisher_OutputFcn, ...
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

function Postpolisher(handles)
%% settings for despiking
    init.despike            =get(handles.despikeon,'Value');            %Despiking on
    init.spikemaxwidth      =str2double(get(handles.width,'String'));  
    init.updownmargin       =str2double(get(handles.margin,'String'));  %fraction that steps can be different
    init.spikeup            =get(handles.dirup,'Value');
    init.spikedown          =get(handles.dirdown,'Value');
    init.spikeboth          =get(handles.dirboth,'Value');
%% settings for slope merging
    init.slopemerge         =get(handles.mergeon,'Value');                  %Merging on
    init.wmin               =str2double(get(handles.widthmerge,'String'));  %Max width
%% Bootstrap settings
    init.booton             =get(handles.erroreston,'Value'); %bootstrapping on
    init.bootoff            =get(handles.errorestoff,'Value'); %bootstrapping off
    Postpolisher_mainloop(init)
    
    function Postpolisher_mainloop(init)
        display("go go go")
        init



% --- Executes just before PostPolisher is made visible.
function PostPolisher_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to PostPolisher (see VARARGIN)

% Choose default command line output for PostPolisher
handles.output = hObject;
pwd=cd;
set(handles.dirboth,'enable','On')
set(handles.dirup,'enable','On')
set(handles.dirdown,'enable','On')
set(handles.width,'enable','On')
set(handles.margin,'enable','On')
set(handles.widthmerge,'enable','Off')
set(handles.errorestoff,'Value',1)
set(handles.erroreston,'Value',0)
set(handles.directory,'String',pwd)
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes PostPolisher wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = PostPolisher_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function widthmerge_Callback(hObject, eventdata, handles)
% hObject    handle to widthmerge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of widthmerge as text
%        str2double(get(hObject,'String')) returns contents of widthmerge as a double


% --- Executes on button press in go.
function go_Callback(~, ~, handles)
% hObject    handle to go (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(PostPolisher);
Postpolisher(handles)


function directory_Callback(hObject, eventdata, handles)
% hObject    handle to directory (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of directory as text
%        str2double(get(hObject,'String')) returns contents of directory as a double



function margin_Callback(hObject, eventdata, handles)
% hObject    handle to margin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of margin as text
%        str2double(get(hObject,'String')) returns contents of margin as a double




function width_Callback(hObject, eventdata, handles)
% hObject    handle to width (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of width as text
%        str2double(get(hObject,'String')) returns contents of width as a double


% --- Executes on button press in despikeon.
function despikeon_Callback(hObject, eventdata, handles)
set(handles.dirboth,'enable','On')
set(handles.dirup,'enable','On')
set(handles.dirdown,'enable','On')
set(handles.width,'enable','On')
set(handles.margin,'enable','On')
set(handles.widthmerge,'enable','Off')

% --- Executes on button press in mergeon.
function mergeon_Callback(hObject, eventdata, handles)
set(handles.dirboth,'enable','Off')
set(handles.dirup,'enable','Off')
set(handles.dirdown,'enable','Off')
set(handles.width,'enable','Off')
set(handles.margin,'enable','Off')
set(handles.widthmerge,'enable','On')