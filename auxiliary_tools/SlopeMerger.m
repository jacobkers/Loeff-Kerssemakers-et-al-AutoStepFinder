function [FitX,step_props]=SlopeMerger(wmin,XX,FitX)
%% Remove slope series by merging
%Aim 1: The Autostepfinder tends to fit trains of small steps to smooth
%slopes of non-instant steps, for example if the data is low-pass filtered.
%This tools finds such step trains and merges
%them by retaining the middle step index. After that, fit is re-done

%Description: 
%1) SlopeMerger: no input, runs a demo trace
%2) SlopeMerger(wmin): load data and process.
    %with this option, output is saved with extension '_merged'
%3) SlopeMerger(wmin,XX,FitX)  

%Output:
 %XX is a single data column, FitX is its fit. 
  %wmin(3) is the maximum distance aqual-signed steps should be spaced to 
  %be considered for merging   
    
%Output: new fit and step properties

%Extras: includes autorun demo option of a  noisy step trace

%Reference: Auxiliairy tool of AutoStepfinder: 
%A fast and automated step detection method for single-molecule analysis.
%Luuk Loeff*, Jacob Kerssemakers*, Chirlmin Joo & Cees Dekker.
% * Equal contribution
% Last update by Jacob Kerssemakers, October 2020
%--------------------------------------------------------------------------

 
 
 
%% Get trace
if (nargin==1)  %Load standard Autostepfinder output
    [FileName,PathName] = uigetfile('*_fits.txt', 'Open the fit of interest');
    data=double(dlmread(strcat(PathName,'/',FileName),',',2,1));
    XX=data(:,1);
    FitX=data(:,2);
end  
 %% Demo: make a simple trace with a number of step trains  or spikes
demo=0;
if nargin==0, demo=1;end
if demo 
    wmin=3;
    [XX,FitX]=MakeSlopedStepsTrace;   
end

% ---------------------------------------------
 close all;
 Lx=length(XX);
 TT=(1:Lx)';
[step_props,~,~]=Get_StepsFromFit(TT,FitX);  
oriFitX=FitX;

%% remove slopes
for rp=1:2
    [~,~,~,together]=find_slopesteps(step_props(:,1), step_props(:,5), wmin);  
    FitX=Get_FitFromSteps(XX,together,'median');
    [step_props,~,~]=Get_StepsFromFit(TT,FitX);
end

if nargin <3
    plot(XX,'k-','LineWidth', 2); hold on;
    plot(oriFitX, 'c-','LineWidth', 1.5);
      plot(FitX, 'r-', 'LineWidth', 2);
     legend('data','original fit','slope-merged fit');
     dum=1;
end

if nargin==1
     Axz=(1:length(XX))';
     data=[Axz XX FitX];
     outname1=strcat(PathName,'/',FileName(1:end-4),'_merged.txt');
     dlmwrite(outname1,data);
     outname2=strcat(PathName,'/',FileName(1:end-4),'_merged_props.txt');
     dlmwrite(outname2,step_props);
end


    


function  [steptrains,isolated,merged,together]=find_slopesteps(index, steps, wmin)
%find all connected trains of steps; output is a structure 'steptrains' and
%sets of indices: 
%'isolated', 
%'merged'  set of one index of thelargest step per step train
%together: 'isolated' and 'merged' combined
%----------------------------------------------


ls=length(steps);

%get indices close to at least one another AND same sign
dt_fw=diff(index);
ds_fw=diff(sign(steps));
sel=find((dt_fw<=wmin)&ds_fw==0);                   %nearby with same sign
if ~isempty(sel)
narrows=unique(sort([index(sel) ; index(sel+1)]));  %these are potential trains
isolated=index(find(~ismember(index,narrows)));     %these are isolated
merged=[];

%proceed with points nearby another

%initialize first set
first_ix2add=narrows(1);                    %first index of set              %
sign2add=sign(steps(index==first_ix2add));  %sign of set
steptrains(1).ix=first_ix2add;  
steptrains(1).sign=sign2add;

set_index=1; 
keep=narrows(find(first_ix2add~=narrows));

keep_steps=steps(find(ismember(index,keep)));    

while length(keep)>0
    set_grows=1;
while set_grows   
    current_ones=steptrains(set_index).ix;
    Lset=length( current_ones);
    new_ones=[];
    for ii=1:Lset  %for all set members
        this_i=steptrains(set_index).ix(ii);
        this_trainsign=steptrains(set_index).sign;
        dt=abs(keep-this_i);                    %time difference
        ds=sign(keep_steps)-this_trainsign;     %sign difference
        sel=find((dt<=wmin)&(ds==0));           %includes same-sign
        %%find one new point per existing to the train ; update remainder:
        if ~isempty(sel)  
            one2add=keep(sel(1));
            new_ones=[new_ones one2add];   
            keep=keep(find(keep~=one2add));  %remove from remainder
            keep_steps=steps(find(ismember(index,keep)));            
        end
    end
    %add the set of new points to the train
    if ~isempty(new_ones)
        new_ones=unique(new_ones);
        steptrains(set_index).ix=unique(sort([current_ones new_ones]));
    else  %set complete; get step closest to middle and go to next set
        this_train=steptrains(set_index).ix';
        these_steps=abs(steps(ismember(this_train,index)));
        %1) find largest step       
        %[~,idx]=max(these_steps);
        %2) closest to middle
        [~,idx]=min(abs(this_train-mean(this_train)));
        
        
        merged_train=this_train(idx);
        merged=[merged ; merged_train];
        %prepare for next
        set_grows=0;
        set_index=set_index+1; 
        %initialize the next step train:
        if ~isempty(keep)
            first_ix2add=keep(1);                 %get one new point
            sign2add=sign(steps(index==first_ix2add));  %sign of set           
            keep=keep(find(keep~=first_ix2add));  %remove from remainder
            keep_steps=steps(find(ismember(index,keep)));
            steptrains(set_index).ix=first_ix2add; %first index of next set
            steptrains(set_index).sign=sign2add;   %sign of set
        end
    end
end
end
together=sort(unique([isolated ; merged]));
else
    steptrains=[];
    merged=[];
    isolated=index;
    together=index;
end

  function [StepsX,levelX, histX]=Get_StepsFromFit(T,FitX)
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
        
    
    function FitX=Get_FitFromSteps(X,indexlist,modus)       
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

function [XX,FitX]=MakeSlopedStepsTrace
 Lxi=1000;
    noize=1;                   %noise level
    XX=[]; FitX=[];
    stepsizes=[-5:2:5];
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
            
            %smooth both
            XX=smooth(XX,10);
            FitX=round(smooth(FitX,10));

         
         