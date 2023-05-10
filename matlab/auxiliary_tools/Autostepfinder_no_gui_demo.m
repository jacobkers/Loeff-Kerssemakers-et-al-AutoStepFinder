function Autostepfinder_no_gui_demo
%% This is a stub that shows how  'AutoStepfinder_no_gui' can be used
%Jacob Kerssemakers 2023.

%get default settings
initval=set_default;

%make a simple example data trace together with its 'ground truth' trace:
SaveName='Demo';
[XX,GroundTruth]=DemoTrace;

%% example of some overrides and call AutoStepfinder:
close all
%the two settings that matter most:
initval.SMaxTreshold    = 0.15; 
initval.fitrange        = 10000;
%control the output:
initval.no_save=1;
initval.no_plot=1;
[FinalFit,FinalSteps]=AutoStepfinder_no_gui(XX, SaveName,initval);
%external demo plot:
if initval.no_plot
    set(gcf, 'units', 'normalized', 'position', [0.02 0.38 0.71 0.5]);
    plot(XX, 'LineWidth',2,'Color',[0,0.2,1]); hold on;
    plot(GroundTruth,   'LineWidth',1.5,  'Color','r');  
    plot(FinalFit,   'LineWidth',1.5,  'Color',[1,0.7,0]);  
    xlabel('Time (s)','FontSize',12);                                   
    ylabel('Position (A.U.)','FontSize',12);                            
    set(gca,'TickDir','out','TickLength',[0.003 0.0035],'box', 'off');
    initval.LegPlt=legend('data','ground truth', 'fit');                       
    set(initval.LegPlt,'box', 'Off','Orientation','Horizontal');
end



function initval=set_default
% Parameters (& default value)
    %%specific for non-gui version to facilitate clean nested usage:
    initval.no_save=0;
    initval.no_plot=0;
    
    %% paths
    initval.datapath        = pwd;          %starting data path    
    initval.codefolder      = pwd;    
    initval.GlobalErrorAccept=0.1;           %User value for accepting a split or merge round solution   
   
    %% essentials
    initval.SMaxTreshold    = 0.15; 
    initval.fitrange        = 1000;       
    
    %% tuning and basic processing
    initval.overshoot       = 1;            %Increase of decrease the number to-be-fitted steps relative to the determined optimum.      
    initval.fitmean         = 1;        %Use mean for fitting, else median
    initval.manualon        = 0;        %Manual mode
    initval.N_manualsteps   = 100;      %if used, this forces a set number of steps   
    
    %% post processing and various
    initval.stepnumber      = initval.fitrange;    %Iteration range of the measurement
    initval.nextfile        = 0;                                
    initval.resolution      = 1;        %Resolution of measurement
    initval.meanbase        = 0;        %Mean value of the base line
    initval.max_range       = 40;       %Max range for noise estimation       
    initval.basetresh     = -1E20;      %base line treshhold level (set -1E20 to turn off)
    
    %% show result for diagnosis
    initval.show_fitplot =      1;
    initval.show_userplot =     1;       
    initval.show_scurve_plot  = 1;      
    
    %% saving formats
    initval.txtoutput       = 1;        %Output .txt files
    initval.matoutput       = 0;        %Output .mat files    
    initval.parametersout   = 1;        %Save parameters output file.
    initval.fitsoutput      = 1;        %Save fits output file.
    initval.propoutput      = 1;        %Save properties output file.
    initval.scurvesoutput   = 1;        %Save S-curves output file.

    %% noise and step error handling
    initval.estimatenoise   = 0;        %data noise estimation on
    initval.bootstraprepeats=0;     % %bootstrap eror cycles per step (0=off)
    
    %% run settings
     %Note: it is advised to leave the 'run' settings as-is as these are
     %linked to the original GUI version. This version is understood to
     %work with one file at a time.        
    initval.singlerun       = 1;        %Single or batch run
    if initval.singlerun    == 1
      initval.hand_load     =  1;       %Single Run
      initval.rerun         =  0;
      if initval.rerun      == 1
         initval.hand_load =  0;        %load last configuration
      end
    else    
      initval.hand_load     =  2;       %Batch Run
      initval.datapath      = uigetdir(initval.datapath);   %Get directory for batch analysis
    end    
    
function [XX,FitX]=DemoTrace
    %Simple Demo Trace
    Lxi=2000;
    noize=5;                   %noise level
    XX=[]; FitX=[];
    stepsizes=[-10 10 30 50 5];
    Nspikes=30; 
    spike_ampli=100;
    isolated=[1 [10:50:Lxi]];  %just a step train
    indices=isolated;
    Li=length(indices);
    for ii=1:Li
        dice=ceil(length(stepsizes)*(rand(1,1)));
        levs(ii)=stepsizes(dice);  %add random step levels
    end
    %build initial trace
    for ii=1:Li-1
        L_sect=indices(ii+1)-indices(ii)+1;
        XX=[XX ; noize*randn(L_sect, 1)+levs(ii)];
        FitX=[FitX ; zeros(L_sect, 1)+levs(ii)];
    end 
