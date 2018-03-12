%This function performs a step-fitting algorithm on data in a quick fashion.
%Background theory can be found in  Supplement 3 of the paper: "Assembly dynamics of microtubules at molecular resolution" 
%J.W. Kerssemakers, E. L. Munteanu, L. Laan, T.L. Noetzel, M.E. Janson, M. Dogterom 
%Nature  442(7103) 709-12 (2006). 
%Algorhitm code developed by Jacob Kerssemakers
%GUI shell developed by Luuk Loeff
%---------------------------------
  
function varargout = Stepfinder_with_gui(varargin)
% STEPFINDER_WITH_GUI MATLAB code for Stepfinder_with_gui.fig
%      STEPFINDER_WITH_GUI, by itself, creates a new STEPFINDER_WITH_GUI or raises the existing
%      singleton*.
%
%      H = STEPFINDER_WITH_GUI returns the handle to a new STEPFINDER_WITH_GUI or the handle to
%      the existing singleton*.
%
%      STEPFINDER_WITH_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in STEPFINDER_WITH_GUI.M with the given input arguments.
%
%      STEPFINDER_WITH_GUI('Property','Value',...) creates a new STEPFINDER_WITH_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the STEPFINDER_WITH_GUI before Stepfinder_with_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Stepfinder_with_gui_OpeningFcn via varargin.
%
%      *See STEPFINDER_WITH_GUI Options on GUIDE's Tools menu.  Choose "STEPFINDER_WITH_GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Stepfinder_with_gui

% Last Modified by GUIDE v2.5 18-Nov-2016 11:57:42

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Stepfinder_with_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @Stepfinder_with_gui_OutputFcn, ...
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

% --- Executes during object creation, after setting all properties.
function data_path_CreateFcn(hObject, eventdata, handles)
% hObject    handle to data_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
data_directory = pwd; %insert your default directory here
%data_directory = 'K:\bn\alg\Shared\Luuk_Jacob\Gui\Testdata';
set(hObject,'String', num2str(data_directory));
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');

end

function StepfinderSuperAuto2016(handles) 
     %% Parameters set in GUI
    initval.datapath        = get(handles.data_path, 'string');         %Data path
    
    initval.codefolder      = pwd;
    
    initval.GlobalErrorAccept=0.1;                                       %User value for accepting a split or merge round solution
    
    initval.SMaxTreshold    = str2double(get(handles.SMaxTreshold,...   %Threshold for second round of fitting
                              'string')); 
    initval.overshoot       = str2double(get(handles.overshoot,...      %Increase of decrease the number to-be-fitted steps relative to the determined optimum.
                              'string'));     
    initval.fitrange        = str2double(get(handles.fitrange,...       %Number of steps to be fitted
                              'string'));       
    initval.stepnumber      = initval.fitrange; 
    initval.nextfile        = 1;
    initval.steprepulsion   = str2double(get(handles.steprepulsion,...  %This term prevents very small steps (2 sample points) to be fitted
                              'string'));                                   
    initval.resolution      = str2double(get(handles.res_mes,...        %Resolution of measurement
                              'string'));      
    initval.meanbase        = str2double(get(handles.meanbase,...       %Mean value of the base line
                              'string'));       
    initval.overbase        = str2double(get(handles.baseover,...       %Baseline overshoot value  
                              'string'));          
    initval.userplt         = get(handles.userplton,'Value');           %Turn user plot function on/ off
    initval.scurve_eval     = get(handles.scurveeval,'Value');          %Turn S-curve evaluation on/ off
    initval.fitmean         = get(handles.fitmean,'Value');             %Use mean for fitting
    initval.fitmedian       = get(handles.fitmedian,'Value');           %Use median for fitting
    initval.treshonoff      = get(handles.basetreshon,'Value');         %Turn base line treshholding on/ off
      if initval.treshonoff == 1
      initval.basetresh     = initval.meanbase*initval.overbase;        %Treshhold the mean of your base line
      else
      initval.basetresh     = -100000;
      end
        
    initval.singlerun       = get(handles.singrun,'Value');             %Single or batch run
      if initval.singlerun  == 1
      initval.hand_load     =  1;                                       %Single Run
      else    
      initval.hand_load     =  2;                                       %Batch Run
      initval.datapath      = uigetdir(initval.datapath);               %Get directory for batch analysis
      end 
           
    %% Remaining parameters --> discuss what to keep.  
    initval.setsteps=0;  %If larger than 0, this will set the numbers of steps to be fitted! 
    initval.CropInputDataFactor=1;
    initval.showintermediateplots=0;    

while initval.nextfile>0;  
    [T,X,SaveName,initval]=Get_Data(initval); % Load data, check for NaN + Inf values    
    stepnumber_firstrun=min([ceil(length(T)/4) initval.fitrange]);        
    ResiduX=X;  FitX=0*X;
    S_Curves=zeros(stepnumber_firstrun+1,2);                   
    AllSteps=[];  
    for fitround=1:2;
        initval.stepnumber=stepnumber_firstrun;                         
        [FitResiduX,~,S_Curve]=StepfinderCore(ResiduX,initval);       
        steproundaccept=(max(S_Curve)>initval.SMaxTreshold);
        if steproundaccept
           [Steps, ~, ~]=Get_StepsFromFit_MeanLevel(T,X,FitResiduX);
            Steps=AddStep_Errors(ResiduX,Steps);  %measured error
            Steps(:,9)=fitround; AllSteps=[AllSteps; Steps];
        end   
        S_Curves(:,fitround)=S_Curve;
        ResiduX=ResiduX-FitResiduX;  %new residu  
        FitX=FitX+FitResiduX;        %new fit
    end  
    if isempty(AllSteps), disp('No steps found'); else  %Final analysis:
      [FinalSteps, FinalFit]=BuildFinalfit_ViaStepErrors(T,X,AllSteps);  

    %% output section
    
    SaveStepsUserFormat(initval,FinalSteps,SaveName);                %Step properties
    %Fits
      Time                      = T*initval.resolution;                 %Time Axis
      Data                      = X;                                    %Raw Data 
      fits_table                = table(Time, Data, FinalFit);          %Save variables in table
      writetable(fits_table, [SaveName,'_fits.txt']);                   %Save table containing fits
    %S-Curve
      Stepnumber                = (1:1:length(S_Curves))';              %Stepnumbers
      SCurveRound1              = S_Curves(:,1);                        %S-Curve round 1
      SCurveRound2              = S_Curves(:,2);                        %S-Curve round 2
      SCurve_table              = table(Stepnumber, SCurveRound1,...    %Save variables in table 
                                  SCurveRound2);        
      writetable(SCurve_table, [SaveName,'_SCurve.txt']);               %Save table containing S-curves     
     
%% Plotting
%Plotting in GUI   
    close(findobj('type','figure','name','S-Curve Evaluation'));        %close S-curve plots --> for batch mode
    close(findobj('type','figure','name','User plots'));                %close user plots --> for batch mode
    cla;                                                                %clear axes 
    axis(handles.plot_fit);
    plot(Time,Data,...                                                  %Plot Data
    'LineWidth',1,....                                                  %Linewidth
    'Color',[0,0.2,1]);                                                 %Color line RBG
    hold on
    plot(Time,FinalFit,...                                              %Plot Fit
    'LineWidth',2,....                                                  %Linewidth
    'Color',[1,0.7,0]);                                                 %Color line RBG
    initval.MaxX=Time(end);                                             %Determine length X axis
    initval.MaxY=max(Data)*1.2;                                         %Determine length Y axis
    initval.MinY=min(Data);                                         %Determine length Y axis
    xlim([0 initval.MaxX]);                                             %Set X axis
    ylim([initval.MinY initval.MaxY]);                                             %Set Y axis
    xlabel('Time (s)','FontSize',12);                                   %Label X axis
    ylabel('Position (A.U.)','FontSize',12);                            %Label Y axis
    set(gca,'TickDir','out','TickLength',[0.003 0.0035],'box', 'off');  %Set ticks outslide plotting area, remove box plot
    initval.LegPlt=legend('Data','Fit');                                        %Set labels legend
    set(initval.LegPlt,'box', 'Off','Orientation','Horizontal');                %Remove box legend, align horizontal
    if initval.treshonoff == 1                                          %If baseline tresholding is on plot line
    initval.BaseLine=repmat(initval.basetresh,1,length(Time));                  %Generate array filled with baseline value
    hold on
    plot(Time,initval.BaseLine,...                                              %Plot treshold line                                              
    'LineWidth',2,...
    'Color',[1,0,0]); 
    initval.LegPlt=legend('Data','Fit','Treshold');                             %Update labels legend
    set(initval.LegPlt,'box', 'Off','Orientation','Horizontal');                %Remove box legend, align horizontal
    end
    
    if initval.singlerun == 0                                           %Save figure with fit, during batch mode
    cd(initval.datapath);                                               %Set current directory to datapath
    initval.FitPlt = figure('visible', 'off');                          %Plot invisible figure for saving
    set(gcf, 'units', 'normalized', 'position', [0.01 1 0.7 0.5])       %Set size figure same as GUI
    copyobj(handles.plot_fit, initval.FitPlt);                          %Copy figure from GUI
    xlabel('Time (s)','FontSize',12);                                   %Set label X axis
    ylabel('Position (A.U.)','FontSize',12);                            %Set label Y axis
    set(gca,'TickDir','out','TickLength',[0.003 0.0035],'box', 'off');  %Set ticks outslide plotting area, remove box plot
    saveas(initval.FitPlt, [SaveName '_Fit.jpg']);                      %Save figure as jpg
    close(initval.FitPlt);                                              %Close invisible plot
    end
    cd(initval.codefolder);                                             %Set current directory to pwd

    
% Plotting user Plots
    if initval.userplt==1
    User_Plot_Result(T,X,FinalFit,FinalSteps,S_Curves,initval);
    if initval.singlerun == 0                                           %Save userplot jpg if batch run is on
    cd(initval.SaveFolder);                                             %Set current directory to savefolder                                                   
    saveas(findobj('type','figure','name','User plots'),...             %Save userplot figure as jpg
    [SaveName '_User_plot.jpg'])
    end
    end

% Plotting S-Curve Evaluation
    if initval.scurve_eval==1                                          
    SCurve_Evaluation(T,X,FinalFit,FinalSteps,S_Curves,initval,ResiduX,FitResiduX,FitX); 
    if initval.singlerun == 0                                           %Save userplot jpg if batch run is on
    cd(initval.SaveFolder);                                             %Set current directory to savefolder  
    saveas(findobj('type','figure','name','S-Curve Evaluation'),...     %Save S-curve figure as jpg
    [SaveName '_SCurve.jpg'])
    end
    end

    disp('steps found:'), display(length(FinalSteps(:,1)))
    end
    disp('done!')
end
     

% --- Executes just before Stepfinder_with_gui is made visible.
function Stepfinder_with_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Stepfinder_with_gui (see VARARGIN)
cla(handles.plot_fit);
axis(handles.plot_fit);
plot(0,0);
set(handles.figure1, 'units', 'normalized', 'position', [0.01 1 0.7 0.5]);
movegui('northwest')
xlabel('Time (s)','FontSize',12);                                       %Do not rotate xlabel
ylabel('Position (A.U.)','FontSize',12, 'rot', 90);                     %Rotate ylabel
set(gca,'TickDir','out','TickLength',[0.003 0.0035],'box', 'off');      %Ticks outslide plotting area
set(handles.PostPros, 'Visible','Off'); 
set(handles.AdvancedSettings,'Visible','Off');
set(handles.AdvancedFitting,'Visible','Off');
set(handles.Scurve_eval,'Visible','Off');
set(handles.advancedoff,'Value',1)
set(handles.Scurve_eval,'Visible','Off')
StepRep = 0.1;
set(handles.steprepulsion, 'String', StepRep)
set(handles.fitmean,'Value',1)
set(handles.scruveevaloff,'Value',1)
set(handles.singrun,'Value',1);
set(handles.userpltoff,'Value',1);

% Choose default command line output for Stepfinder_with_gui
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Stepfinder_with_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Stepfinder_with_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit1_Callback(hObject, eventdata, ~)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function data_path_Callback(hObject, eventdata, handles)
% hObject    handle to data_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of data_path as text
%        str2double(get(hObject,'String')) returns contents of data_path as a double


function fitrange_Callback(hObject, eventdata, handles)
% hObject    handle to fitrange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fitrange as text
%        str2double(get(hObject,'String')) returns contents of fitrange as a double


% --- Executes during object creation, after setting all properties.
function fitrange_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fitrange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function res_mes_Callback(hObject, eventdata, handles)
% hObject    handle to res_mes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of res_mes as text
%        str2double(get(hObject,'String')) returns contents of res_mes as a double


% --- Executes during object creation, after setting all properties.
function res_mes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to res_mes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function SMaxTreshold_Callback(hObject, eventdata, handles)
% hObject    handle to SMaxTreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SMaxTreshold as text
%        str2double(get(hObject,'String')) returns contents of SMaxTreshold as a double


% --- Executes during object creation, after setting all properties.
function SMaxTreshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SMaxTreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function steprepulsion_Callback(hObject, eventdata, handles)
% hObject    handle to steprepulsion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of steprepulsion as text
%        str2double(get(hObject,'String')) returns contents of steprepulsion as a double


% --- Executes during object creation, after setting all properties.
function steprepulsion_CreateFcn(hObject, eventdata, handles)
% hObject    handle to steprepulsion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function overshoot_Callback(hObject, eventdata, handles)
% hObject    handle to overshoot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of overshoot as text
%        str2double(get(hObject,'String')) returns contents of overshoot as a double


% --- Executes during object creation, after setting all properties.
function overshoot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to overshoot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in runprogram.
function runprogram_Callback(hObject, eventdata, handles,Time,Data,FinalFit)
% hObject    handle to runprogram (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
StepfinderSuperAuto2016(handles)


% --- Executes on button press in treshholdbox.
function treshholdbox_Callback(hObject, eventdata, handles)
% hObject    handle to treshholdbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of treshholdbox



function meanbase_Callback(hObject, eventdata, handles)
% hObject    handle to meanbase (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of meanbase as text
%        str2double(get(hObject,'String')) returns contents of meanbase as a double


% --- Executes during object creation, after setting all properties.
function meanbase_CreateFcn(hObject, eventdata, handles)
% hObject    handle to meanbase (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on key press with focus on treshholdbox and none of its controls.
function treshholdbox_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to treshholdbox (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


function baseover_Callback(hObject, eventdata, handles)
% hObject    handle to baseover (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of baseover as text
%        str2double(get(hObject,'String')) returns contents of baseover as a double


% --- Executes during object creation, after setting all properties.
function baseover_CreateFcn(hObject, eventdata, handles)
% hObject    handle to baseover (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when selected object is changed in uibuttongroup1.
function uibuttongroup1_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uibuttongroup1 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
initval.PostProcessOn=get(handles.basetreshon,'Value');
if initval.PostProcessOn == 1
       set(handles.PostPros, 'Visible','On');
end
initval.PostProcessOff=get(handles.basetreshoff,'Value');
if initval.PostProcessOff == 1
       set(handles.PostPros, 'Visible','Off');   
       MeanBase = 0;
       set(handles.meanbase, 'String', MeanBase);
       BaseOver = 1.10;
       set(handles.baseover, 'String', BaseOver);
end


% --- Executes when selected object is changed in paneladv.
function paneladv_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in paneladv 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
initval.AdvancedOn=get(handles.advancedon,'Value');
if initval.AdvancedOn == 1
       set(handles.AdvancedSettings,'Visible','On')
       set(handles.AdvancedFitting,'Visible','On')
       set(handles.Scurve_eval,'Visible','On')
end
initval.AdvancedOff=get(handles.advancedoff,'Value');
if initval.AdvancedOff == 1
       set(handles.AdvancedSettings,'Visible','Off')
       set(handles.AdvancedFitting,'Visible','Off')
       set(handles.Scurve_eval,'Visible','Off')
       StepRep = 0.1;
       set(handles.steprepulsion, 'String', StepRep)
       set(handles.fitmean,'Value',1)
       set(handles.scruveevaloff,'Value',1)
end

% --- Executes when figure1 is resized.
function figure1_SizeChangedFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



























%% This section (570 to ~770) contains the 'Core' function of the stepfinder; 
%it can be cut and autorun independently (on a simple simulated curve) for demo purposes
function [FitX,stepsX,S_fin]=StepfinderCore(X,initval)
%This function splits data in a quick fashion.
%This one is a compact version of the first quick 2007 version
%output: list of stepsizes: [index time  levelbefore levelafter step dwelltimeafter steperror]
if nargin<2
    X=2*round(0.5+0.1*sin(2*pi*(1:20000)'/500))+rand(20000,1); 
    initval.stepnumber=300; initval.overshoot=1;   
end
 tic 
 disp('number of data points'); disp(length(X))
    %% 1 split, estimate best fit, repeat
    initval.stepnumber=min([ceil(length(X)/4) initval.stepnumber]);
    [~,~,S_raw]=Split_until_ready(X,initval); %run 1: full iteration
    
    [bestshot,S_fin]=Eval_Scurve(S_raw);

    initval.stepnumber=round(min([( bestshot-1)*initval.overshoot ceil(length(X)/4)])); 
    [FitX,~,~]=Split_until_ready(X,initval); %run2: estimate done by the program
    stepsX=Get_Steps(FitX); 

toc   
if nargin<2
    close all;
    subplot(2,1,1); plot(X,'r'),hold;plot(FitX,'k','LineWidth',2);
    title('Data and Fit');xlabel('time');ylabel('position,a.u.');
    subplot(2,1,2); semilogx(S_fin,'-o');
    title('S-curve');xlabel('Stepnumber');ylabel('S-value, a.u.');
end


function [bestshot,S_fin]=Eval_Scurve(S_raw);
    S_raw(S_raw<1)=1; S2=S_raw-1;  %remove base line
    BaseLine=linspace(0,S2(end),length(S2)); S3=S2-BaseLine';
    %S4=smooth(S3,ceil(ix/25));
    [~,i1]=max(S3);  %peak
    %sel=find(S3>0.9*pk1);i2=max(sel(sel>=i1)); bestshot=i2;
    bestshot=i1;
    S_fin=S3;
  
    
function stepsX=Get_Steps(FitX)
%list of stepsizes: [index time step levelbefore levelafter dwelltimeafter]    
    lx=length(FitX);
    T=(1:lx)';
    difX=FitX(2:lx)-FitX(1:lx-1);
    sel=find(difX~=0);    
    lsel=length(sel);
    dwellX=T(sel(2:lsel))-T(sel(1:lsel-1)); dwellX=[dwellX' T(lx)-T(sel(lsel))]';
    stepsX=[sel T(sel) FitX(sel) FitX(sel+1) difX(sel) dwellX]; 
            
function [FitX,f,S]=Split_until_ready(X,initval)
     c=1; stop=0;
     N=length(X);    
     FitX=mean(X)*ones(N,1); 
     S=ones(initval.stepnumber,1);

     %Create the first plateau------------------------------------------
     istart=1; istop=length(X);
     [inxt, avl, avr,rankit]=Splitfast(X(istart:istop));           
     f=[[1, 1, 1, 0, 0,0];
        [istart, istop, inxt+istart-1, avl, avr,rankit]; ...
        [N, N, N,0, 0,0]];  
     %parameters needed for calculating S(1):-----------------
    qx=sum(X.^2);                                   %sum of squared data
    qm=N*(mean(X))^2;                               %sum of squared averages plateaus, startvalue
    aqm=(inxt-istart+1)*avl^2+(istop-inxt)*avr^2;   %sum of squared averages anti-plateaus, startvalue
    S(c)=(qx-aqm)/(qx-qm);                          %S: ratio of variances of fit and anti-fit        
    %---------------------------------       

     while stop==0; %Split until ready----------------------------------
        c=c+1;
        fsel=find((f(:,2)-f(:,1)>5)&f(:,6)~=0);        %among those plateaus sensibly long..
        [~,idx2]=max(f(fsel,6)); idx=(fsel(idx2));   %...find the best candidate to split. 
        FitX=Adapt_Fit(f,idx,FitX);                     %adapt fit-curve
        [f,qm,aqm]=expand_f(f,qm,aqm,idx,X);            %adapt plateau-table; adapt S
        S(c)=(qx-aqm)/(qx-qm);                             %Calculate new S-function  
        stop=(1.0*c>initval.stepnumber);
    end   %-------------------------------------------------------------------
          
function [f,qm,aqm]=expand_f(f,qm,aqm,idx,X)
%this function inserts two new plateau-property rows on the location of one old one
%....and adapts the S-function nominator and denominator

   %1) Label new levels.FLR locatess 'plateau fit  left right' etc  
    nFLR=f(idx-1,2)-f(idx-1,3);       avFLR=f(idx-1,5);      %FLR
    nFML=f(idx,3)-f(idx,1)+1;         avFML=f(idx,4);        %FML
    nFMR=f(idx,2)-f(idx,3);           avFMR=f(idx,5);        %FMR
    nFRL=f(idx+1,3)-f(idx+1,1)+1;     avFRL=f(idx+1,4);      %FRL
      
    %remove contribution from old plateau(s) from S-function terms
    qm=qm-1/(nFML+nFMR)*(nFML*avFML+nFMR*avFMR)^2;            %FM
    aqm=aqm-    1/(nFLR+nFML)*(nFLR*avFLR+nFML*avFML)^2-...   %CL
                1/(nFMR+nFRL)*(nFMR*avFMR+nFRL*avFRL)^2;      %CR
            
    %2a) construct new first plateau entry, left
	istart=f(idx,1); istop=f(idx,3);
	[inxt, avl, avr,rankit]=Splitfast(X(istart:istop));
    n1=[istart istop inxt+istart-1, avl, avr,rankit];
    
    %2b) construct new first plateau entry, right
	istart=f(idx,3)+1; istop=f(idx,2);
	[inxt, avl, avr,rankit]=Splitfast(X(istart:istop));
	n2=[istart istop inxt+istart-1, avl, avr,rankit];
	
    %3) Insert these two new plateaus in place of the old one
	[lf,~]=size(f);      
    block1=f(1:idx,:);      block1(idx,:)=n1;
	block2=f(idx:lf,:); 	block2(1,:)=n2;
	f=[block1; block2];
    
    %Label newly defined levels.FMLR locates 'plateaufit mid/left/right' etc
    nFMLL=f(idx,3)-f(idx,1)+1;          avFMLL=f(idx,4);    %FMLL
    nFMLR=f(idx,2)-f(idx,3);            avFMLR=f(idx,5);    %FMLR
    nFMRL=f(idx+1,3)-f(idx+1,1)+1;      avFMRL=f(idx+1,4);  %FMRL
    nFMRR=f(idx+1,2)-f(idx+1,3);        avFMRR=f(idx+1,5);  %FMRR
    
    %4) add contribution from new plateau(s) to S-function terms
    qm=qm...
        +1/(nFMLL+nFMLR)*(nFMLL*avFMLL+nFMLR*avFMLR)^2 ...  %FML
        +1/(nFMRL+nFMRR)*(nFMRL*avFMRL+nFMRR*avFMRR)^2;     %FMR
    aqm=aqm ...
        +1/(nFLR+nFMLL)*(nFLR*avFLR+nFMLL*avFMLL)^2 ...      %CTL
        +1/(nFMLR+nFMRL)*(nFMLR*avFMLR+nFMRL*avFMRL)^2 ...   %CTM
        +1/(nFMRR+nFRL)*(nFMRR*avFMRR+nFRL*avFRL)^2;         %CTR

function FitX=Adapt_Fit(f,idx,FitX)
	%This function creates step  and property curves adds new plateaus
	i1=f(idx,1); i2=f(idx,3);av1=f(idx,4);
    i3=f(idx,3)+1; i4=f(idx,2);av2=f(idx,5);
    FitX(i1:i2)=av1; FitX(i3:i4)=av2;
              
 function [idx, avl, avr,rankit]=Splitfast(Segment)              %
%this function also adresses a one-dim array 'Segment'
%and determines the best step-fit there
%To save time, functions like 'mean' are avoided
    w=length(Segment);     
    if w>3
		Chisq=(1:w-1)*0;                           
        AvL=Segment(1);    AvR=sum(Segment(2:w))/(w-1); AvAll=sum(Segment)/w;  
        for t=2:w-2
            AvL=(AvL*(t-1)+Segment(t))/t;     AvR=(AvR*(w-t+1)-Segment(t))/(w-t);
            DStepL=AvL-AvAll;           DStepR=AvR-AvAll;
            DeltaChisq=((DStepL.^2)*t+(DStepR.^2)*(w-t));            
            Chisq(t)=-DeltaChisq;       
        end
         [~,idx]=min(Chisq(2:w-2)); idx=idx+1;
         avl=mean(Segment(1:idx));                    avr=mean(Segment(idx+1:w));
         %rankit=(avr-avl)^2/(1/(idx-1)+1/(w-idx-1));  %quantity expresing predicted accuracy step (squared)
         rankit=(avr-avl)^2*w;
    else                                            %low-limit cases
         rankit=0;               
         switch w
             case 3 
                a=(Segment(2)+Segment(3))/2-Segment(1);            b=Segment(3)-(Segment(1)+Segment(2))/2;
                cL=[Segment(1) , (Segment(1)+Segment(2))/2];      cR=[(Segment(2)+Segment(3))/2 , Segment(3)]; 
                [~,idx]=max([a b]);    avl=cL(idx);  avr=cR(idx);
            case 2
                idx=1;  avl=Segment(1); avr=Segment(2); 
        end
    end


    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
%% This section contains code related to the multipass steps
function StepsX=AddStep_Errors(X,StepsX)   
%This function calculates the errors associated with the steps via the
%standard deviation of the surrounding plateaus
[ls,col]=size(StepsX); i1=0;
for i=1:ls
    i2=StepsX(i);
    if i<ls
        i3=StepsX(i+1);
    else 
        i3=length(X);
    end
    rmsbefore=std(X(i1+1:i2)); 
    Nbefore=i2-i1;
    rmsafter=std(X(i2+1:i3)) ;
    Nafter=i3-i2;
    StepsX(i,col+1)=2*(rmsbefore^2/Nbefore+rmsafter^2/Nafter)^0.5; %plus minus 95%
    i1=i2;
end

function [FinalSteps, FinalFit]=BuildFinalfit_ViaStepErrors(T,X,AllSteps)
%build a step fit based on all retained indices. Perform step-by-step error
%analysis to accept second-round (residual) steps or not (first-round steps are always
%accepted at this stage)
%1) first, add 
    [~,ix]=sort(AllSteps(:,1));  %sort by index (ignoring round)
    AllSteps=AllSteps(ix,:);
    RoundNo=AllSteps(:,9);
    CandidateFit=Get_FitFromStepsindices(X,AllSteps(:,1),'mean');
    CandidateSteps=Get_StepsFromFit_MeanLevel(T,X,CandidateFit); 
    LC=length(CandidateSteps(:,1));
    CandidateRoundNo=zeros(LC,1);
    for ii=1:LC
        idxC=CandidateSteps(ii,1);
        sel=find(AllSteps(:,1)==idxC);
        CandidateRoundNo(ii)=RoundNo(sel(1));
    end
    CandidateSteps=AddStep_Errors(X,CandidateSteps);
    CandidateRelStepErr=(CandidateSteps(:,8)./abs(CandidateSteps(:,5))); 
    [~,~,FinalErrorTreshold]=Outlier_flag(CandidateRelStepErr,2,0.8,'positive',0);
    

%% keep round 1-steps AND 'good' round 2 steps
    sel=find((CandidateRelStepErr<2*FinalErrorTreshold)|CandidateRoundNo==1);  
    FinalIdxes=CandidateSteps(sel,1);
    FinalFit=Get_FitFromStepsindices(X,FinalIdxes,'mean');
    FinalSteps=Get_StepsFromFit_MeanLevel(T,X,FinalFit); 
    FinalSteps=AddStep_Errors(X,FinalSteps);
    LF=length(FinalSteps(:,1));
    FinalRoundNo=zeros(LF,1);
    for ii=1:LF
        idxC=FinalSteps(ii,1);
        sel=find(CandidateSteps(:,1)==idxC);
        FinalRoundNo(ii)=CandidateRoundNo(sel(1));
    end        
    FinalSteps(:,9)=FinalRoundNo;   
   
 
function [T,X,SaveName,initval]=Get_Data(initval);

    %This function loads the data, either standard or user-choice
    disp('Loading..');
    CurrentFolder=pwd;
    switch initval.hand_load     
        case 1      
        cd(initval.datapath);
        [FileName,PathName] = uigetfile('*.*','Select the signal file');
        cd(CurrentFolder);
        source=strcat(PathName,FileName);
        data=double(dlmread(source));
        SaveName=FileName(1:length(FileName)-4);
        initval.nextfile=0;      
        case 2
            fileindex=initval.nextfile;
            cd(initval.datapath);
            AllFileNames=dir('*.txt');
            cd(CurrentFolder);
            AllFiles=length(AllFileNames);
            FileName=AllFileNames(fileindex).name;
            data=double(dlmread(strcat(initval.datapath,'\',FileName)));
            SaveName=FileName(1:length(FileName)-4);
            if initval.nextfile==AllFiles
                initval.nextfile=0;
            else
                initval.nextfile=fileindex+1;
            end
            
    end
    LD=length(data(:,1));
    if initval.CropInputDataFactor<1
        Cropit=round(LD/initval.CropInputDataFactor);    
        T=data(1:Cropit,1); 
        X=data(1:Cropit,2); 
    else
        T=data(:,1); 
        X=data(:,2); 
    end
%    check for NaN + Inf values
     infcheck=isinf(X);
     nancheck=isnan(X);

     if sum(infcheck)>0;
         msgbox('Your data contains infinite values','ERROR', 'error')
         return;
     end

     if sum(nancheck)>0;
         msgbox('Your data contains NaN values','ERROR', 'error')
         return;
     end      
     T=(1:1:length(X))';   
     disp('Analyzing:'), disp(SaveName);

  function FitX=Get_FitFromStepsindices(X,indexlist,modus) 
  % This function builds plateau data
    %list of levels: [startindex  stopindex starttime stoptime level dwell stepbefore stepafter]
        lx=length(X);
        lsel=length(indexlist); %note: index points to last point before step

        %Build a 'FitX' based on the median levels (ot on the averages)
        idxes=[0 ; indexlist ; lx];
        FitX=0*X;
        for ii=1:lsel+1
            ixlo=idxes(ii)+1;  %first index of plateau
            ixhi=idxes(ii+1);  %last index
            switch modus
                case 'mean', FitX(ixlo:ixhi)=nanmean(X(ixlo:ixhi));
                case 'median', FitX(ixlo:ixhi)=nanmedian(X(ixlo:ixhi));
            end
        end    
        
        
 function [StepsX,levelX, histX]=Get_StepsFromFit_MeanLevel(T,X,FitX);
%This function builds tables of steps or levels properties from a step fit
%Values are based on Averages of plateaus
    lx=length(FitX);
    difX=FitX(2:lx)-FitX(1:lx-1);
    sel=find(difX~=0);  %note: index points to last point before step
    lsel=length(sel);
        
    %dwell time after
    dwellT=T(sel(2:lsel))-T(sel(1:lsel-1)); 
    dwellTafter=[dwellT' T(lx)-T(sel(lsel))]';
    dwellTbefore=[T(sel(1))-T(1) dwellT']'; 
    StepsX=[sel T(sel) FitX(sel) FitX(sel+1) difX(sel) dwellTbefore dwellTafter]; 
    %list of stepsizes: [index time levelbefore levelafter step dwelltime before dwelltimeafter]
    
    Levstartidx=[1; sel+1];  %first index of level
    Levstopidx=[sel; lx];    %last index of level
    LevstartT=T(Levstartidx);
    LevstopT=T(Levstopidx);
    LevLevel=[FitX(sel); FitX(lx)];
    LevDwell=Levstopidx-Levstartidx+1;
    LevStepBefore=[0; difX(sel)];
    LevStepAfter=[difX(sel); 0];
    levelX=[Levstartidx Levstopidx LevstartT LevstopT LevLevel LevDwell LevStepBefore LevStepAfter];
    %list of levels: [startindex starttime stopindex stoptime level dwell stepbefore stepafter]
    
    
    %histogram
    stepsizes=StepsX(:,4)-StepsX(:,3);
    mx=max(stepsizes); mn=min(stepsizes);
    stp=(mx-mn)/20;
    hx=(mn:stp:mx);
    if ~isempty(hx)
        histX=hist(stepsizes,hx)';
        histX=[hx' histX];
    else
        histX=[1 1];
    end
    
    
 function [MedFitX, MedStepsX, MedLevelX, MedHistX]=Get_StepsFromFit_MedianLevel(T,X,FitX) 
%This function calculates plateau medians  based on step locations;  
%Steps: list of stepsizes: [index time levelbefore levelafter step dwelltime before dwelltimeafter]
%Levels: [startindex  stopindex starttime stoptime level dwell stepbefore stepafter]
    
    lx=length(FitX);
    difX=FitX(2:lx)-FitX(1:lx-1);
    sel=find(difX~=0);  %note: index points to last point before step
    lsel=length(sel);
 
    %Build a 'FitX' based on the median levels (ot on the averages)
    idxes=[0 ; sel ; lx];
    MedFitX=0*X;
    for ii=1:lsel+1
        ixlo=floor(idxes(ii)+1);  %first index of plateau
        ixhi=floor(idxes(ii+1));  %last index
        MedFitX(ixlo:ixhi)=nanmedian(X(ixlo:ixhi));
    end
       
    difMedX=MedFitX(2:lx)-MedFitX(1:lx-1);
    
    %dwell time after
    dwellT=T(sel(2:lsel))-T(sel(1:lsel-1)); 
    dwellTafter=[dwellT' T(lx)-T(sel(lsel))]';
    dwellTbefore=[T(sel(1))-T(1) dwellT']'; 
    MedStepsX=[sel T(sel) MedFitX(sel) MedFitX(sel+1) difMedX(sel) dwellTbefore dwellTafter]; 
    %list of stepsizes: [index time levelbefore levelafter step dwelltime before dwelltimeafter]
    
    Levstartidx=[1; sel+1];  %first index of level
    Levstopidx=[sel; lx];    %last index of level
    LevstartT=T(Levstartidx);
    LevstopT=T(Levstopidx);
    LevLevel=[MedFitX(sel); MedFitX(lx)];
    LevDwell=Levstopidx-Levstartidx+1;
    LevStepBefore=[0; difMedX(sel)];
    LevStepAfter=[difMedX(sel); 0];
    MedLevelX=[Levstartidx Levstopidx LevstartT LevstopT LevLevel LevDwell LevStepBefore LevStepAfter];
    %list of levels: [startindex  stopindex starttime stoptime level dwell stepbefore stepafter]
    
    %histogram
    stepsizes=MedStepsX(:,4)-MedStepsX(:,3);
    mx=max(stepsizes); mn=min(stepsizes);
    stp=(mx-mn)/50;
    hx=(mn:stp:mx);
    if ~isempty(hx)
        MedHistX=hist(stepsizes,hx)';
        MedHistX=[hx' MedHistX];
    else
        MedHistX=[1 1];
    end
    
    function [flag,cleandata,treshold]=Outlier_flag(data,tolerance,sigchange,how,sho);
%this function is meant to find a representative value for a standard
%deviation in a heavily skewed distribution (typically, flat data with
% %peaks). It calculates the standard deviation and average the data;
% Based on these, outliers are determined and excluded for a new calculation
% of average and SD; this is repeated until sigma does not change anymore too much
% . This is repeated until the new sigma does not change much
% %anymore
%output: positions of outliers

%Jacob Kers 2013 and before---------------------------------------------
binz=50;


if nargin<5  %For testing/demo purposes
    close all
    data=JK00_DEMODATA_Peaks;
    tolerance=2;
    sigchange=0.7;
    how='positive';
    sho=1;
    plot(data,'o-');
    binz=20;
end

sigma=1E20;            %at start, use a total-upper-limit 
ratio=0;
ld=length(data);
flag=ones(ld,1);  %at start, all points are selected
cleandata=data;
while ratio<sigchange     %if not too much changes anymore; the higher this number the less outliers are peeled off.
    sigma_old=sigma;
    selc=find(flag==1);
    data(flag==1); 
    ls=length(selc);
    av=nanmedian(data(selc));       %since we expect skewed distribution, we use the median iso the mea     
    sigma=nanstd(data(selc));
    ratio=sigma/sigma_old;
    treshold=tolerance*sigma+av;
    switch how
        case 'positive',  flag=(data-av)<tolerance*sigma;     %adjust outlier flags
        case 'all',  flag=abs(data-av)<tolerance*sigma;     %adjust outlier flags  
    end
    %plot menu------------------  
    if sho==1
        cleandata=data(selc); 
        hx=(min(cleandata):(range(cleandata))/binz:max(cleandata));   %make an axis
        sthst=hist(cleandata,hx);
        figure;
        bar(hx,sthst);
        title('Histogram');
        dum=ginput(1);
        pause(0.5);  
        close(gcf);
    end
    %---------------------------- 
end
cleandata=data(selc); 
hx=(min(cleandata):(range(cleandata))/binz:max(cleandata));   %make an axis
sthst=hist(cleandata,hx);

function SaveStepsUserFormat(initval,FinalSteps,SaveName)
%This function saves data in a user-specific way.
 disp('Saving Files...')

     initval.SaveFolder                 =  [SaveName,'_Fitting_Result'];        %Make new folder to save results
     initval.SaveFolder                 =  fullfile(initval.datapath, initval.SaveFolder);
         if ~exist(initval.SaveFolder, 'dir')                                   %Check if folder already exists
         mkdir(initval.SaveFolder);                                             %If not create new folder
         end
         cd(initval.SaveFolder);                                                %Set current directory to new folder
     %Step properties
      idx_steps                 =  FinalSteps(:,4)> initval.basetresh;  %Treshholding of steps
      IndexStep                 =  FinalSteps(idx_steps,1);             %Index where step occured
      TimeStep                  =  FinalSteps(idx_steps,2)...           %Time when step occured
                                   *initval.resolution; 
      LevelBefore               =  FinalSteps(idx_steps,3);             %Level before step
      LevelAfter                =  FinalSteps(idx_steps,4);             %Level after step
      StepSize                  =  LevelAfter - LevelBefore;            %Size of step
      DwellTimeStepBefore       =  FinalSteps(idx_steps,6)...           %Dwelltime step before
                                   *initval.resolution;       
      DwellTimeStepAfter        =  FinalSteps(idx_steps,7)...           %Dwelltime step after
                                   *initval.resolution;  
      StepError                 =  FinalSteps(idx_steps,8);             %Error of each step
      
      properties_table          = table(IndexStep,TimeStep,...          %Save variables in table
                                  LevelBefore,LevelAfter,StepSize,...
                                  DwellTimeStepBefore,DwellTimeStepAfter,StepError); 
                              
      writetable(properties_table, [SaveName,'_properties.txt']);       %Save table containing properties
 
 function SCurve_Evaluation(T,X,FinalFit,FinalSteps,S_Curves,initval,Xpeel,FitXpeel, FitX); 
    Stepnumber     = (1:1:length(S_Curves))'; 
    StepRange=max(Stepnumber);
    idx_steps= FinalSteps(:,4)> initval.basetresh; % Treshholding of steps.
    Time           = T*initval.resolution; %Time Axis
    SCurveRound1   = S_Curves(:,1); %S-Curve round 1
    SCurveRound2   = S_Curves(:,2); %S-Curve round 2

figure('Name','S-Curve Evaluation','NumberTitle','off','units', 'normalized', 'position', [0.745 0.32 0.25 0.6]);
%S-Curve round 1
    SCurve1Plt=subplot(3,1,1);
    cla(SCurve1Plt);
    SmoothFact = (initval.fitrange/1000);
    SmoothCurve1 = smooth(SCurveRound1,SmoothFact);
    IdxMaxiCurve1 = find(max(SCurveRound1) == SCurveRound1);
    xMaxCurve1 = Stepnumber(IdxMaxiCurve1);
    yMaxCurve1 = SmoothCurve1(IdxMaxiCurve1);
    plot(Stepnumber,SmoothCurve1,'k-');
    xlim([0 initval.fitrange]); 
    set(gca,'TickDir','out','TickLength',[0.003 0.0035],'box', 'off'); %Ticks outslide plotting area
    hold on;
    scatter(xMaxCurve1,yMaxCurve1,'bo')
    title('Round 1');
    xlabel('Step Number');
    ylabel('S-Score');
    
%S-Curve round 2    
    SCurve2Plt=subplot(3,1,2);
    cla(SCurve2Plt);
    SmoothFact = (initval.fitrange/1000);
    IdxMaxiCurve2 = find(max(SCurveRound2) == SCurveRound2);
    SmoothCurve2 = smooth(SCurveRound2, SmoothFact);
    xMaxCurve2 = Stepnumber(IdxMaxiCurve2);
    yMaxCurve2 = SmoothCurve2(IdxMaxiCurve2);
    yMinCurve2 = min(IdxMaxiCurve2);
    plot(Stepnumber,SmoothCurve2,'k-');
    set(gca,'TickDir','out','TickLength',[0.003 0.0035],'box', 'off'); %Ticks outslide plotting area
    hold on;
    scatter(xMaxCurve2,yMaxCurve2,'bo')
    title('Round 2');
    xlabel('Step Number');
    ylabel('S-Score');
    xlim([0 initval.fitrange]);
    ylimCurve2 = yMaxCurve2;
    ylim([yMinCurve2 ylimCurve2]);
    hold on;
    if initval.SMaxTreshold > (yMaxCurve2)
        xtext=(max(initval.fitrange)/2);
        ytext=yMaxCurve2*1.4;
        textRound2='Accuracy was set above the maximum';
        %textRound2=[textRound2 newline 'Second fitting round was not peformed'];
        %text(xtext,ytext,textRound2, 'HorizontalAlignment','center', 'Color', [1,0,0]);
    else
    TreshHold=repmat(initval.SMaxTreshold,1,length(Stepnumber));
    plot(Stepnumber,TreshHold,...
    'LineWidth',1,...
    'Color',[1,0,0]);
    end    
%S-Curve round 2  
    FitPeel=subplot(3,1,3);
    cla(FitPeel)
    plot(Time, Xpeel,...
    'LineWidth',1,....
    'Color',[0,0.2,1]); 
    hold on;
    plot(Time, FitXpeel,...
    'LineWidth',2,....
    'Color',[1,0.7,0]); 
    hold on;
    MaxX=Time(end);
    xlim([0 MaxX]);
    xlabel('Time (s)','FontSize',12);
    ylabel('Position (A.U.)','FontSize',12);
    set(gca,'TickDir','out','TickLength',[0.003 0.0035],'box', 'off'); %Ticks outslide plotting area

     
      
      
  function User_Plot_Result(T,X,FinalFit,FinalSteps,S_Curves,initval); 
%This function can be used to present user (and experiment-specific plots
%using the standard stepfinder output
          
Stepnumber     = (1:1:length(S_Curves))'; 
StepRange=max(Stepnumber);
idx_steps= FinalSteps(:,4)> initval.basetresh; % Treshholding of steps.
SCurveRound1   = S_Curves(:,1); %S-Curve round 1
SCurveRound2   = S_Curves(:,2); %S-Curve round 2
IndexStep      =  FinalSteps(idx_steps,1); %Index where step occured
TimeStep       =  FinalSteps(idx_steps,2)*initval.resolution; %Time when step occured
LevelBefore    =  FinalSteps(idx_steps,3); %Level before step
LevelAfter     =  FinalSteps(idx_steps,4); %Level after step
StepSize       =  LevelAfter - LevelBefore; %Size of step
DwellTimeStep  =  FinalSteps(idx_steps,7)*initval.resolution; %Dwelltime step  
    

figure('Name','User plots','NumberTitle','off','units', 'normalized', 'position', [0.01 0.05 0.3 0.3]);
% %%  Transition Density Plot
%     DensityPlt = subplot(2,1,1);
%     cla(DensityPlt);
%     DensityBinsize = 0.02;
%     Xaxis=0:DensityBinsize:1;
%     Yaxis=0:DensityBinsize:1;
%     BinC={(Xaxis+(0.5*DensityBinsize)),((Yaxis+0.5*DensityBinsize))};
%     values = hist3([LevelAfter(:), LevelBefore(:)],'Ctrs',BinC);
%     contour(Xaxis,Yaxis,values);
%     colorbar;
%     title('Transition Density Plot');
%     xlabel('Level Before');
%     ylabel('Level After');

%Step-Level
    StepLPlt = subplot(2,1,1);
    cla(StepLPlt);
    StepBinsize = 0.05;
    XStep=0:StepBinsize:1;
    [histcounts, ~]= hist(LevelAfter,XStep);
    bar(XStep,histcounts)
    title('Levels');
    xlabel('Step Level');
    ylabel('Counts');
    ylimstep=max(histcounts*1.1);
    ylim([0 ylimstep]);
    xlim([0 1]);
        
%Step-size
    StepSPlt = subplot(2,1,2);
    cla(StepSPlt);
    StepBinsize = 0.05;
    XStep=-1:StepBinsize:1;
    [histcounts, ~]= hist(StepSize,XStep);
    bar(XStep,histcounts)
    title('Step-size');
    xlabel('Step Size');
    ylabel('Counts');
    ylimstep=max(histcounts*1.1);
    ylim([0 ylimstep]);
    xlim([-1 1]);
    
    
 
     








