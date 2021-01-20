
function varargout = StepMaker(varargin);
% 

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @StepMaker_OpeningFcn, ...
                   'gui_OutputFcn',  @StepMaker_OutputFcn, ...
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


function Stepmaker(handles)
%GUI parameters
initval.flatstep        =get(handles.flatstep,'Value');              %Flat distr steps
initval.minstep         =str2double(get(handles.minstep,'String'));
initval.maxstep         =str2double(get(handles.maxstep,'String'));
initval.gausstep        =get(handles.gausstep,'Value');              %Gaus distr steps
initval.meanstep        =str2double(get(handles.meanstep,'String'));
initval.sigmastep       =str2double(get(handles.sigmastep,'String'));
initval.expstep         =get(handles.expstep,'Value');               %Exp distr steps
initval.decaystep       =str2double(get(handles.decaystep,'String'));
initval.flatdwell       =get(handles.flatdwell,'Value');             %Flat dwell dwell
initval.mindwell        =str2double(get(handles.mindwell,'String'));
initval.gausdwell       =get(handles.gausdwell,'Value');             %Gausd distr dwell
initval.maxdwell        =str2double(get(handles.maxdwell,'String'));
initval.meandwell       =str2double(get(handles.meandwell,'String'));
initval.expdwell        =get(handles.expdwell,'Value');              %Flat distr dwell
initval.sigmadwell      =str2double(get(handles.sigmadwell,'String'));
initval.decaydwell      =str2double(get(handles.decaydwell,'String'));
initval.stepsnumber     =str2double(get(handles.stepsnumber,'String'));
initval.noisesteps      =str2double(get(handles.noisesteps,'String'));
initval.repeats         =str2double(get(handles.repeats,'String'));
initval.addbase         =str2double(get(handles.addbase,'String'));
initval.repeatsteps     =str2double(get(handles.traces,'String'));

Stepmakermainloop(initval,handles)

    
function Stepmakermainloop(initval,handles)
        display("Lets make steps!!!!")



% --- Executes just before StepMaker is made visible.
function StepMaker_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to StepMaker (see VARARGIN)
set(handles.flatstep,'Value',1);              %Flat distr steps
set(handles.minstep,'String', -10);
set(handles.maxstep,'String',10);
set(handles.gausstep,'Value',0);              %Gaus distr steps
set(handles.meanstep,'String', 10);
set(handles.sigmastep,'String', 3);
set(handles.expstep,'Value', 0);               %Exp distr steps
set(handles.decaystep,'String',100);
set(handles.flatdwell,'Value',1);             %Flat dwell dwell
set(handles.mindwell,'String', -10);
set(handles.gausdwell,'Value', 0);             %Gausd distr dwell
set(handles.maxdwell,'String',100);
set(handles.meandwell,'String',75);
set(handles.expdwell,'Value',0);              %Flat distr dwell
set(handles.sigmadwell,'String', 25);
set(handles.decaydwell,'String', 100);
set(handles.stepsnumber,'String', 20);
set(handles.noisesteps,'String', 3);
set(handles.repeats,'String', 1);
set(handles.addbase,'String', 0);
set(handles.traces,'String', 5);
set(handles.maxstep, 'Enable','On');
set(handles.minstep, 'Enable','On');
set(handles.meanstep, 'Enable','Off');
set(handles.sigmastep, 'Enable','Off');
set(handles.decaystep, 'Enable','Off');
set(handles.maxdwell, 'Enable','On');
set(handles.mindwell, 'Enable','On');
set(handles.meandwell, 'Enable','Off');
set(handles.sigmadwell, 'Enable','Off');
set(handles.decaydwell, 'Enable','Off');
% Choose default command line output for StepMaker
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);





% --- Outputs from this function are returned to the command line.
function varargout = StepMaker_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function stepsnumber_Callback(hObject, ~, ~)
checkmax_stepsnumber=get(hObject,'String');
checkmax_stepsnumber=isnan(str2double(checkmax_stepsnumber));
     if checkmax_stepsnumber==1
         msgbox('The number of steps is NaN.','ERROR', 'error')
         set(hObject,'String',20);
     return;
     end
if checkmax_stepsnumber < 1
         msgbox('The number of steps is smaller than 1. The input value has been set to 1','ERROR', 'error')
         set(hObject,'String',1);
     return;     
end 
 

function noisesteps_Callback(hObject, ~, ~)
checkmax_noisesteps=get(hObject,'String');
checkmax_noisesteps=isnan(str2double(checkmax_noisesteps));
     if checkmax_noisesteps==1
         msgbox('The noise setting is NaN.','ERROR', 'error')
         set(hObject,'String',3);
     return;
     end

   
function minstep_Callback(hObject, ~, ~)
checkmax_minstep=get(hObject,'String');
checkmax_minstep=isnan(str2double(checkmax_minstep));
     if checkmax_minstep==1
         msgbox('The Min stepsize setting is NaN.','ERROR', 'error')
         set(hObject,'String',-10);
     return;
     end


function GenerateData_Callback(hObject, ~, ~)
% hObject    handle to GenerateData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Stepmaker(handles)

function decaystep_Callback(hObject, ~, ~)
checkmax_decaystep=get(hObject,'String');
checkmax_decaystep=isnan(str2double(checkmax_decaystep));
     if checkmax_decaystep==1
         msgbox('The decay setting is NaN.','ERROR', 'error')
         set(hObject,'String',100);
     return;
     end



function mindwell_Callback(hObject, ~, ~)
checkmax_mindwell=get(hObject,'String');
checkmax_mindwell=isnan(str2double(checkmax_mindwell));
     if checkmax_mindwell==1
         msgbox('The Min dwell time setting is NaN.','ERROR', 'error')
         set(hObject,'String',50);
     return;
     end
          if checkmax_mindwell < 1
         msgbox('The min dwell time is smaller than 1. The input value has been set to 1','ERROR', 'error')
         set(hObject,'String',1);
     return;     
end 
   



function meandwell_Callback(hObject, ~, ~)
checkmax_meandwell=get(hObject,'String');
checkmax_meandwell=isnan(str2double(checkmax_meandwell));
     if checkmax_meandwell==1
         msgbox('The mean dwell time setting is NaN.','ERROR', 'error')
         set(hObject,'String',75);
     return;
     end
if checkmax_meandwell < 2
         msgbox('The number of steps is smaller than 1. The input value has been set to 1','ERROR', 'error')
         set(hObject,'String',2);
     return;     
end 



function decaydwell_Callback(hObject, ~, ~)
checkmax_decaydwell=get(hObject,'String');
checkmax_decaydwell=isnan(str2double(checkmax_decaydwell));
     if checkmax_decaydwell==1
         msgbox('The decay setting is NaN.','ERROR', 'error')
         set(hObject,'String',100);
     return;
     end




function maxdwell_Callback(hObject, ~, ~)
checkmax_maxdwell=get(hObject,'String');
checkmax_maxdwell=isnan(str2double(checkmax_maxdwell));
     if checkmax_maxdwell==1
         msgbox('The max dwell time setting is NaN.','ERROR', 'error')
         set(hObject,'String',100);
     return;
     end
     if checkmax_maxdwell < 2
         msgbox('The max dwelltime is smaller than 2. The input value has been set to 2','ERROR', 'error')
         set(hObject,'String',2);
     return;     
end 





function sigmadwell_Callback(hObject, ~, ~)
checkmax_sigmadwell=get(hObject,'String');
checkmax_sigmadwell=isnan(str2double(checkmax_sigmadwell));
     if checkmax_sigmadwell==1
         msgbox('The sigma dwell time setting is NaN.','ERROR', 'error')
         set(hObject,'String',25);
     return;
     end



function addbase_Callback(hObject, ~, ~)
checkmax_addbase=get(hObject,'String');
checkmax_addbase=isnan(str2double(checkmax_addbase));
     if checkmax_addbase==1
         msgbox('The add baseline setting is NaN.','ERROR', 'error')
         set(hObject,'String',0);
     return;
     end

function meanstep_Callback(hObject, ~, ~)
checkmax_meanstep=get(hObject,'String');
checkmax_meanstep=isnan(str2double(checkmax_meanstep));
     if checkmax_meanstep==1
         msgbox('The mean stepsize setting is NaN.','ERROR', 'error')
         set(hObject,'String',10);
     return;
     end

function sigmastep_Callback(hObject, ~, ~)
checkmax_sigmastep=get(hObject,'String');
checkmax_sigmastep=isnan(str2double(checkmax_sigmastep));
     if checkmax_sigmastep==1
         msgbox('The sigma stepsize setting is NaN.','ERROR', 'error')
         set(hObject,'String',10);
     return;
     end



function repeats_Callback(hObject, ~, ~)
checkmax_repeats=get(hObject,'String');
checkmax_repeats=isnan(str2double(checkmax_repeats));
     if checkmax_repeats==1
         msgbox('The number of repeats is NaN.','ERROR', 'error')
         set(hObject,'String',1);
     return;
     end
     if checkmax_repeats < 1
         msgbox('The number of repeats is smaller than 1. The input value has been set to 1','ERROR', 'error')
         set(hObject,'String',1);
     return;     
end 
 


function traces_Callback(hObject, ~, ~)
checkmax_traces=get(hObject,'String');
checkmax_traces=isnan(str2double(checkmax_traces));
     if checkmax_traces==1
         msgbox('The number of traces is NaN.','ERROR', 'error')
         set(hObject,'String',1);
     return;
     end
if checkmax_traces < 1
         msgbox('The number of traces is smaller than 1. The input value has been set to 1','ERROR', 'error')
         set(hObject,'String',1);
     return;     
end 
 

function maxstep_Callback(hObject, ~, ~)
checkmax_maxstep=get(hObject,'String');
checkmax_maxstep=isnan(str2double(checkmax_maxstep));
     if checkmax_maxstep==1
         msgbox('The max stepsize setting is NaN.','ERROR', 'error')
         set(hObject,'String',1);
     return;
     end


% --- Executes on button press in flatstep.
function flatstep_Callback(~, ~, handles)
% hObject    handle to flatstep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.maxstep, 'Enable','On');
set(handles.minstep, 'Enable','On');
set(handles.meanstep, 'Enable','Off');
set(handles.sigmastep, 'Enable','Off');
set(handles.decaystep, 'Enable','Off');



% --- Executes on button press in gausstep.
function gausstep_Callback(~, ~, handles)
set(handles.maxstep, 'Enable','Off');
set(handles.minstep, 'Enable','Off');
set(handles.meanstep, 'Enable','On');
set(handles.sigmastep, 'Enable','On');
set(handles.decaystep, 'Enable','Off');



% --- Executes on button press in expstep.
function expstep_Callback(~, ~, handles)
set(handles.maxstep, 'Enable','Off');
set(handles.minstep, 'Enable','Off');
set(handles.meanstep, 'Enable','Off');
set(handles.sigmastep, 'Enable','Off');
set(handles.decaystep, 'Enable','On');


% --- Executes on button press in flatdwell.
function flatdwell_Callback(~, ~, handles)
set(handles.maxdwell, 'Enable','On');
set(handles.mindwell, 'Enable','On');
set(handles.meandwell, 'Enable','Off');
set(handles.sigmadwell, 'Enable','Off');
set(handles.decaydwell, 'Enable','Off');


% --- Executes on button press in gausdwell.
function gausdwell_Callback(~, ~, handles)
set(handles.maxdwell, 'Enable','Off');
set(handles.mindwell, 'Enable','Off');
set(handles.meandwell, 'Enable','On');
set(handles.sigmadwell, 'Enable','On');
set(handles.decaydwell, 'Enable','Off');


% --- Executes on button press in expdwell.
function expdwell_Callback(~, ~, handles)
set(handles.maxdwell, 'Enable','Off');
set(handles.mindwell, 'Enable','Off');
set(handles.meandwell, 'Enable','Off');
set(handles.sigmadwell, 'Enable','Off');
set(handles.decaydwell, 'Enable','On');
