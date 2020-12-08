function StepMaker(tracetype)
%JWJK_B:-------------------------------------------------------------------
%Title:Simulate Step Traces
%Summary: %this function  creates a simple stepped trace with noise. 
%Approach: Steps  are taken as random samples from a 
%stepsize - and dwell time distribution. 
%This distribution can be refined in the 'Set_type_of_signal' function

%Dwell times can be set as 'dependent' on the step size or 'independent'.
%In the latter case they are taken as random samples from a flat or
%exponential dwell time distribution.

%By varying the settings, a wide range of traces types can be generated

%Below, example settings are given for case studies

%Input: none
%Output: step trace saved as text column in root path

%Jacob Kerssemakers, Cees Dekker Lab, 2018-20
%JWJK_B:-------------------------------------------------------------------

%% Case examples
%See line 223 and below how these work
if nargin==0
    all_examples=[
     {'DownTrendGaussian'};
     {'DownStepsWithSpikes'};
     {'FlatDistribution'};
     {'StepTrainDown'};%
     {'BleachTrace'};
     {'Custom User'}]     
else
     all_examples=[{tracetype}];
end
for ii=1:length(all_examples)
    example_signal=char(all_examples{ii});
    init=Set_type_of_signal(example_signal);
    %:JWJK_B-------------------------------------------------------------------         
        save_it=1;   
        writepath=[pwd, '\test_traces'];
        if ~isdir(writepath), mkdir(writepath);,end
      if strcmp( example_signal,'BleachTrace')
           [data,SaveName]=MakeBleachTrace(init);
      else    
           SaveName=strcat('testdata_', example_signal,...
            init.TypeOfSteps,num2str(init.NumberOfDwells-1),...
            'Smz_',num2str(init.SmoothWindow),'.txt');
        %% distribution generation
        DwellTimes=zeros(init.NumberOfDwells+1,1); 
        binz=50;
        if strcmp(init.StepDwellRelation, 'InDependent')       
            MinDwellTime=1;
            MaxDwellTime=300;

            DwellBinAxis=linspace(MinDwellTime,MaxDwellTime,binz);
            switch init.TypeOfDwells
                case 'Flat'
                    DwellTimeCurve=(1+0*(DwellBinAxis)); 
                case 'Exponential'
                    DwellTimeCurve=exp(-DwellBinAxis/(0.5*MaxDwellTime)); 
            end
            DwellTimeCurve=DwellTimeCurve/sum(DwellTimeCurve);  %normalized dstribution
            SumDwellCurve=cumsum(DwellTimeCurve);
        end

        %% Choose step distributions
        [StepSizeCurve,StepBinAxis,StepDwells]=Set_distributions(init.TypeOfSteps,binz,init); 
        StepSizeCurve=StepSizeCurve/sum(StepSizeCurve);  %normalized dstribution
        SumStepSizeCurve=cumsum(StepSizeCurve); 

      %% Build lists. done by picking a random sample from the step and dwell distributions
        % Build step size list  
        StepSizes=zeros(init.NumberOfDwells,1);  
        for jj=1:init.NumberOfDwells  
                Sample=Pick_from_distribution(StepBinAxis,SumStepSizeCurve);
                StepSizes(jj)=Sample;
        end 

       % Build dwell-time list   
        switch init.StepDwellRelation
            case 'Dependent'  %pick related dwell time
                for jj=1:init.NumberOfDwells
                    st=StepSizes(jj);
                        [~,idx]=min(abs(StepBinAxis-st));
                        DwellTimes(jj)=StepDwells(idx);
                end                  
                case 'InDependent'
                    for jj=1:init.NumberOfDwells
                        Sample=Pick_from_distribution(DwellBinAxis,SumDwellCurve);
                        DwellTimes(jj)=round(Sample);
                    end
       end

        %% Build curve     
        Curve=[];
        Level=0;
        for jj=1:init.NumberOfDwells
            Dwell=DwellTimes(jj);
            Step=StepSizes(jj);        
            NewSection=Level+Step*ones(Dwell,1);        
            NewSection=add_spikes(NewSection,Dwell,init.spike_density,StepSizes,jj,init.NumberOfDwells);
            Curve=[Curve ; NewSection];
            Level=Curve(end);
        end  

        %% Add start or tail
        if init.add_start_or_stop_level<0
            Curve=[zeros(abs(init.add_start_or_stop_level),1)+Curve(1) ; Curve];
        end
        if init.add_start_or_stop_level>0
            Curve=[Curve ; zeros(abs(init.add_start_or_stop_level),1)+Curve(end)];
        end
        Ld=length(Curve);
        TimAx=(1:Ld)';


        %% Add Noise
        MeanStep=median(abs(StepSizes));
        Noise=randn(Ld,1)*init.SignalNoiseRatio*MeanStep;
        Trace=Curve+Noise;
        if init.SmoothWindow>0
            Trace=smooth(Trace,round(init.SmoothWindow));
        end    
        data=[TimAx Trace];

      end   

        %% plot menu
        if nargin<1
            close all;
            figure(298);
            if ~strcmp( example_signal,'BleachTrace')
                subplot(2,2,1);
                bar(StepBinAxis,StepSizeCurve);
                title('Step Sizes');
                subplot(2,2,3); 
                if strcmp(init.StepDwellRelation, 'InDependent')          
                    bar(DwellTimeCurve,'r');
                    title('Dwell Times-independent');
                    xlabel('dwell time');
                    ylabel('occurence, a.u.');
                end
                if strcmp(init.StepDwellRelation, 'Dependent')          
                    bar(StepBinAxis,StepDwells,'r');
                    title('Dwell Times-dependent');
                    xlabel('step size');
                    ylabel('dwell time');
                end
                subplot(1,2,2);
            end

                plot(data(:,1),data(:,2));
                title('Trace')
                legend(strcat(init.TypeOfSteps, '-Distribution'));
                pause(1);
        end

        %% save menu
        if save_it  
            dlmwrite(strcat(writepath,'/',SaveName), data);
        end
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
            spike_duration=3; 
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
            
            
  
   function [StepSizeCurve,StepBinAxis,StepDwells]=Set_distributions(TypeOfSteps,binz,init)
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
            stepsizes=init.Specials.DiscreteSteps.stepsizes;  %sizes
            stepcount=init.Specials.DiscreteSteps.stepcount;    %occurence
            dwelltimes=init.Specials.DiscreteSteps.dwelltimes; %duration
            StepSizeCurve=(0*(StepBinAxis));
            StepDwells=round(0*StepBinAxis);
            for jj=1:length(stepsizes)
                [~,idx]=min(abs(stepsizes(jj)-StepBinAxis));
                 StepSizeCurve(idx)=stepcount(jj);
                 StepDwells(idx)=dwelltimes(jj);
            end
        case 'Gaussian'
            peakpos=init.Specials.forGauss.peakpos;  %in units of the step axis; 0.5 means mid-range
            peakwidth=init.Specials.forGauss.peakwidth; %same
            Rng=MaxStepSize-MinStepSize;
            StepSizeCurve=exp(-( (StepBinAxis-(peakpos*Rng+MinStepSize))/(peakwidth*Rng)).^2);
            StepDwells=round(0*StepBinAxis+10+100*rand(1,binz));
    end
    
    
   function init=Set_type_of_signal(example_signal)
    %% general default step properties 
    %...and explanation
    %(specific settings per type of curve), to be adapted below
    init.TypeOfSteps='DiscreteSteps';  
    %choices are: 'DiscreteSteps' ,  'FlatDistribution', 'Gaussian'
    
    init.NumberOfDwells=20;
    %Note: this is actually the number of plateaus; steps is one lower
    
    init.spike_density=3;  
    %in average events per dwell time; set to 0 to inactivate
    
    init.add_start_or_stop_level=100;  
    %number of time points to add to front or back of curve;
    %negative: to front. %positive: to end
         
    init.StepDwellRelation='Dependent'; 
    %choices are:  'Dependent'; 'InDependent'
    %if Dependent, dwell times are a function of stepsize (as specified in
    %the function 'get_distributions'. If Independent, dwelltimes are 
    %independently treated as a distribution (defined below)
    
    init.TypeOfDwells='Exponential';  
    %choices are: 'Exponential' 'Flat'
    %only applied when using 'Independent'
    
    init.SignalNoiseRatio=(10)^-1;  
    %ratio of noise amplitude relative to average absolute step size
    
    init.SmoothWindow=0;
    %add some low-pass smoothing for checking effects on fit;
    %0=no smoothing
 
    %% if an existing example case was defined, certain settings are adapted:
    switch example_signal
        case 'FlatDistribution'
        % 1) flat distribution:
        % This creates randomly picked steps going up and down
            init.TypeOfSteps='FlatDistribution'; 
            init.NumberOfSteps=10;
            init.spike_density=0; 
            init.add_start_or_stop_level=0; 
            init.StepDwellRelation='Dependent'; 
            init.SignalNoiseRatio=(10)^-1;  
        case 'StepTrainDown'            
            %2) fixed-size steps going down
            init.TypeOfSteps='DiscreteSteps'; 
            init.NumberOfSteps=25;
            init.spike_density=0; 
            init.add_start_or_stop_level=0; 
            init.StepDwellRelation='Dependent'; 
            init.SignalNoiseRatio=(10)^-1;  
            init.Specials.DiscreteSteps.stepsizes=[-10];  %sizes
            init.Specials.DiscreteSteps.stepcount=[1];    %occurence
            init.Specials.DiscreteSteps.dwelltimes=[100]; %duration
            %  This creates a fixed-size step train going down with an end tail
            %-------------------------------------------------------------------------
        case 'DownStepsWithSpikes'         
            % 3) fixed-size steps going down, varying duration, with spikes
            init.TypeOfSteps='DiscreteSteps'; 
            init.NumberOfSteps=10;
            init.spike_density=0.8; 
            init.add_start_or_stop_level=1000; 
            init.StepDwellRelation='InDependent';
            init.TypeOfDwells='Exponential';
            init.SignalNoiseRatio=(10)^-1;  
            init.Specials.DiscreteSteps.stepsizes=[-10];  %sizes
            init.Specials.DiscreteSteps.stepcount=[1];    %occurence
            init.Specials.DiscreteSteps.dwelltimes=[100]; %duration
            %  This creates a step train of varying duration going down with an end
            %  tail, with some blinking events
        case 'DownTrendGaussian'           
            % gaussian distribution, trending down
            init.TypeOfSteps='Gaussian'; 
            init.NumberOfSteps=50;
            init.spike_density=0; 
            init.add_start_or_stop_level=-1000; 
            init.StepDwellRelation='InDependent';
            init.TypeOfDwells='Exponential';
            init.SignalNoiseRatio=(10)^-1;  
            init.Specials.forGauss.peakpos=0.4; %in units of the step axis; 0.5 means mid-range
            init.Specials.forGauss.peakwidth=0.3; %in units of the step axis; peakpos=0.3;
        case 'Custom User'           
            % gaussian distribution, trending down
            init.TypeOfSteps='DiscreteSteps' %'Gaussian';  %'DiscreteSteps'
            init.NumberOfSteps=50;
            init.spike_density=0; 
            init.add_start_or_stop_level=1000; 
            init.StepDwellRelation='Dependent'; %'InDependent'; %'Dependent'; 
            init.TypeOfDwells='Flat';  %'Exponential' 'Flat'
            init.SignalNoiseRatio=(10)^-1;
            init.SmoothWindow=0;
            
            init.Specials.forGauss.peakpos=0.4; %in units of the step axis; 0.5 means mid-range
            init.Specials.forGauss.peakwidth=0.3; %in units of the step axis; peakpos=0.3;
            init.Specials.DiscreteSteps.stepsizes= [-10 -5   5 10];  %sizes
            init.Specials.DiscreteSteps.stepcount= [1    5  1   1 ];    %occurence
            init.Specials.DiscreteSteps.dwelltimes=[200  50 50 200]; %duration
        case 'BleachTrace'
            init.pts=1000; 
            init.tau=init.pts/10; 
            init.N0=7; 
            init.noise=0.3;
    end
    
    function [data,SaveName]=MakeBleachTrace(init);
            %Make a simple bleach-like trace
            tt=linspace(1,1000, init.pts);
            yy=round(init.N0*exp(-tt/init.tau))+init.noise*randn(1,init.pts);
            data=[tt', yy'];
            SaveName=strcat('testdata_simulated_BleachTrace',num2str(init.N0),'steps.txt');
            

  