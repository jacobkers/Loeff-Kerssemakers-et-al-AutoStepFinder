function StepMaker
%JWJK_B:-------------------------------------------------------------------
%Title:Simulate Step Traces
%Summary: %this function  creates a simple stepped trace with noise. 
%Approach: Steps  are taken as random samples from a 
%stepsize - and dwell time distribution. 
%This distribution can be refined in the function 'Set_distributions'
%Dwell times can be set as 'dependent' on the step size or 'independent'.
%In the latter case they are taken as random samples from a flat or
%exponential dwell time distribution.
%By varying the settings, a wide range of traces types can be generated

%Below, example settings are given for case studies

%Input: none
%Output: step trace saved as text column in root path

%Jacob Kerssemakers, Cees Dekker Lab, 2018-20

%% Case examples
%--------------------------------------------------------------------------

%% 1) flat distribution:
%     TypeOfSteps='FlatDistribution'; 
%     NumberOfSteps=10;
%     spike_density=0; 
%     add_start_or_stop_level=0; 
%     StepDwellRelation='Dependent'; 
%     SignalNoiseRatio=(10)^-1;  
%     This creates randomly picked steps going up and down
%-------------------------------------------------------------------------

%% 2) fixed-size steps going down
%     TypeOfSteps='DiscreteSteps'; 
%     NumberOfSteps=10;
%     spike_density=0; 
%     add_start_or_stop_level=1000; 
%     StepDwellRelation='Dependent'; 
%     SignalNoiseRatio=(10)^-1;  
%  ...in additon, go to the function 'Set_distributions' . 
%     Here, under case 'DiscreteSteps'', set:
%             stepsizes=[-10];  %sizes
%             stepcount=[1];    %occurence
%             dwelltimes=[100]; %duration
%  This creates a fixed-size step train going down with an end tail
%-------------------------------------------------------------------------

%% 3) fixed-size steps going down, varying duration, with spikes
%     TypeOfSteps='DiscreteSteps'; 
%     NumberOfSteps=10;
%     spike_density=0.8; 
%     add_start_or_stop_level=1000; 
%     StepDwellRelation='InDependent';
%     TypeOfDwells='Exponential';
%     SignalNoiseRatio=(10)^-1;  
%  ...in additon, go to the function 'Set_distributions' . 
%     Here, under case 'DiscreteSteps'', set:
%             stepsizes=[-10];  %sizes
%             stepcount=[1];    %occurence
%             dwelltimes=[100]; %duration
%  This creates a step train of varying duration going down with an end
%  tail, with some blinking events
%-------------------------------------------------------------------------
%% 4) gaussian distribution, trending down
%     TypeOfSteps='Gaussian'; 
%     NumberOfSteps=50;
%     spike_density=0; 
%     add_start_or_stop_level=-1000; 
%     StepDwellRelation='InDependent';
%     TypeOfDwells='Exponential';
%     SignalNoiseRatio=(10)^-1;  
%  ...in additon, go to the function 'Set_distributions' . 
%     Here, under case 'Gaussian'', set:
%             peakpos=0.3;  %in units of the step axis; 0.5 means mid-range
%             peakwidth=0.5 %in units of the step axis;
%  


%:JWJK_B-------------------------------------------------------------------         
    save_it=1;   
    writepath=pwd;
    
    %% general step properties (specific settings per type of curve)
    TypeOfSteps='DiscreteSteps';  
    %choices are: 'DiscreteSteps' ,  'FlatDistribution', 'Gaussian'
    
    NumberOfDwells=20;
    %Note: this is actually the number of plateaus; steps is one lower
    
    spike_density=0;  
    %in average events per dwell time; set to 0 to inactivate
    
    add_start_or_stop_level=100;  
    %number of time points to add to front or back of curve;
    %negative: to front. %positive: to end
         
    StepDwellRelation='Dependent'; 
    %choices are:  'Dependent'; 'InDependent'
    %if Dependent, dwell times are a function of stepsize (as specified in
    %the function 'get_distributions'. If Independent, dwelltimes are 
    %independently treated as a distribution (defined below)
    
    TypeOfDwells='Exponential';  
    %choices are: 'Exponential' 'Flat'
    %only applied when using 'Independent'
    
    SignalNoiseRatio=(10)^-1;  
    %ratio of noise amplitude relative to average absolute step size
  
    
    %% distribution generation
    DwellTimes=zeros(NumberOfDwells+1,1); 
    binz=50;
    if strcmp(StepDwellRelation, 'InDependent')       
        MinDwellTime=1;
        MaxDwellTime=300;
        
        DwellBinAxis=linspace(MinDwellTime,MaxDwellTime,binz);
        switch TypeOfDwells
            case 'Flat'
                DwellTimeCurve=(1+0*(DwellBinAxis)); 
            case 'Exponential'
                DwellTimeCurve=exp(-DwellBinAxis/(0.5*MaxDwellTime)); 
        end
        DwellTimeCurve=DwellTimeCurve/sum(DwellTimeCurve);  %normalized dstribution
        SumDwellCurve=cumsum(DwellTimeCurve);
    end
   
    %% Choose step distributions
    [StepSizeCurve,StepBinAxis,StepDwells]=Set_distributions(TypeOfSteps,binz); 
    StepSizeCurve=StepSizeCurve/sum(StepSizeCurve);  %normalized dstribution
    SumStepSizeCurve=cumsum(StepSizeCurve); 
    
  %% Build lists. done by picking a random sample from the step and dwell distributions
    % Build step size list  
    StepSizes=zeros(NumberOfDwells,1);  
    for jj=1:NumberOfDwells  
            Sample=Pick_from_distribution(StepBinAxis,SumStepSizeCurve);
            StepSizes(jj)=Sample;
    end 
    
   % Build dwell-time list   
    switch StepDwellRelation
        case 'Dependent'  %pick related dwell time
            for jj=1:NumberOfDwells
                st=StepSizes(jj);
                    [~,idx]=min(abs(StepBinAxis-st));
                    DwellTimes(jj)=StepDwells(idx);
            end                  
            case 'InDependent'
                for jj=1:NumberOfDwells
                    Sample=Pick_from_distribution(DwellBinAxis,SumDwellCurve);
                    DwellTimes(jj)=round(Sample);
                end
   end
    
    %% Build curve     
    Curve=[];
    Level=0;
    for jj=1:NumberOfDwells
        Dwell=DwellTimes(jj);
        Step=StepSizes(jj);        
        NewSection=Level+Step*ones(Dwell,1);        
        NewSection=add_spikes(NewSection,Dwell,spike_density,StepSizes,jj,NumberOfDwells);
        Curve=[Curve ; NewSection];
        Level=Curve(end);
    end  
    
    %% Add start or tail
    if add_start_or_stop_level<0
        Curve=[zeros(abs(add_start_or_stop_level),1)+Curve(1) ; Curve];
    end
    if add_start_or_stop_level>0
        Curve=[Curve ; zeros(abs(add_start_or_stop_level),1)+Curve(end)];
    end
    Ld=length(Curve);
    TimAx=(1:Ld)';
    
    
    %% Add Noise
    MeanStep=median(abs(StepSizes));
    Noise=randn(Ld,1)*SignalNoiseRatio*MeanStep;
    data=[TimAx Curve+Noise];
    
    
    
    %% plot menu
    if nargin<1
        %close all;
        figure(298);
        subplot(2,2,1);
        bar(StepBinAxis,StepSizeCurve);
        title('Step Sizes');
        subplot(2,2,3); 
        if strcmp(StepDwellRelation, 'InDependent')          
            bar(DwellTimeCurve,'r');
            title('Dwell Times-independent');
            xlabel('dwell time');
            ylabel('occurence, a.u.');
        end
        if strcmp(StepDwellRelation, 'Dependent')          
            bar(StepBinAxis,StepDwells,'r');
            title('Dwell Times-dependent');
            xlabel('step size');
            ylabel('dwell time');
            
        end
        subplot(1,2,2);
            plot(data(:,1),data(:,2));
            title('Trace')
            legend(strcat(TypeOfSteps, '-Distribution'));
    end
    
    %% save menu
    if save_it
    SaveName=strcat('testdata_simulated',TypeOfSteps,num2str(NumberOfDwells-1),'.txt');
    dlmwrite(strcat(writepath,'/',SaveName), data);
    end
    
 function Sample=Pick_from_distribution(BinAxis,SumCurve)
            %pick from distribution
            pickval=rand(1);  %pick a random percentile
            sel_st=find(SumCurve>0);
            NonZeroSumCurve=SumCurve(sel_st);
            NonZeroBinAxis=BinAxis(sel_st);                        
            %nearest occupied bin            
            [~,idx]=min(abs(NonZeroSumCurve-pickval));
            %will return first index of plateau which is correct
            Sample=NonZeroBinAxis(idx);  
    
    function NewSection=add_spikes(NewSection,Dwell,spike_density,StepSizes,jj,NumberOfSteps)
        %optional spike adding
        if spike_density>0&(jj>1)&(jj<NumberOfSteps-1)
            spike_duration=5; 
            overlength=1000*Dwell;
            no2add=round(spike_density*1000); %should be larger than 1
            id0=sort(ceil((overlength-1)*rand(no2add,1)));  %two indices)
            sel=find(id0<Dwell-2*spike_duration-1);
            if ~isempty(sel)
                spikestarts=id0(sel); spikestops=id0(sel)+spike_duration;
                for sp=1:length(spikestarts)  %add next step-spike
                    NewSection(spikestarts(sp):spikestops(sp))=...
                    NewSection(spikestarts(sp):spikestops(sp))+StepSizes(jj+1);
                end
                dum=1;
            end
        end         
            
            
  
   function [StepSizeCurve,StepBinAxis,StepDwells]=Set_distributions(TypeOfSteps,binz)
    %set up global range of stepsizes. If 'dependent' dwell times are used,
    %the dwell time is a function of hte stepsize and is specified here. If
    %dwell times are independent, this fucntional dependendence is not used
    %and dwell times are picked from their own distribution.    
    %1) set bin axis 
    MinStepSize=-10;
    MaxStepSize=10; 
    StepBinAxis=linspace(MinStepSize,MaxStepSize,binz);
   
    %2) choose step types
    switch TypeOfSteps
       case 'FlatDistribution'
           %each step in the specified range has an equal chance of
           %occurring
            StepSizeCurve=(1+0*(StepBinAxis));   
            StepDwells=round(0*StepBinAxis+100);
       case 'DiscreteSteps'
            %A limited number of stepsizes is possible; 
            %The 'counts' specify a relative occurence.
            %using only one step gives a step train
            stepsizes=[-10 10];  %sizes
            stepcount=[1 1];    %occurence
            dwelltimes=[50 10]; %duration
            StepSizeCurve=(0*(StepBinAxis));
            StepDwells=round(0*StepBinAxis);
            for jj=1:length(stepsizes)
                [~,idx]=min(abs(stepsizes(jj)-StepBinAxis));
                 StepSizeCurve(idx)=stepcount(jj);
                 StepDwells(idx)=dwelltimes(jj);
            end
        case 'Gaussian'
            peakpos=0.3;  %in units of the step axis; 0.5 means mid-range
            peakwidth=0.5 %same
            Rng=MaxStepSize-MinStepSize;
            StepSizeCurve=exp(-( (StepBinAxis-(peakpos*Rng+MinStepSize))/(peakwidth*Rng)).^2);
            StepDwells=round(0*StepBinAxis+10+100*rand(1,binz));
    end