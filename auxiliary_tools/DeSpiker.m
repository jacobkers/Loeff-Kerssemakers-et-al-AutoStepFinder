function [FitX_merged,StepFitProps_out]=DeSpiker(spikemaxwidth, updownmargin, spikesign,XX,FitX)
%% Remove spikes
%Aim: spikes (short up-and-down) or blinks (reverse) may be of less
%interest to a user. This tool finds and removes these. 

%Description: 
%1) DeSpiker: no input, runs a demo trace on trains
%4) DeSpiker(spikemaxwidth, updownmargin, spikesign): load data and process according to modus.
    %with this option, output is saved with extension '_merged'
%3) DeSpiker(spikemaxwidth, updownmargin, spikesign,XX,FitX)

%Input:
    %XX is a single data column, FitX is its fit.
    %spikemaxwidth(5), : max width for spike to consider
    %updownmargin(0.30): relative magnitude margin per up and down step 
    %spikesign(0): 
%         1: remove positive spikes, 
%        -1: remove negative spikes or 'blinks'
%         0 remove all signs


%Output: new fit and step properties

%Extras: includes autorun demo option of a  noisy step trace

%Reference: Auxiliairy tool of AutoStepfinder: 
%A fast and automated step detection method for single-molecule analysis.
%Luuk Loeff*, Jacob Kerssemakers*, Chirlmin Joo & Cees Dekker.
% * Equal contribution
% Last update by Jacob Kerssemakers, October 2020
%--------------------------------------------------------------------------


%% Get trace
if (nargin==3)  %Load standard Autostepfinder output
    [FileName,PathName] = uigetfile('*_fits.txt', 'Open the fit of interest');
    data=double(dlmread(strcat(PathName,'/',FileName),',',2,1));
    XX=data(:,1);
    FitX=data(:,2);
end  
 %% Demo: make a simple trace with a number of  spikes
if nargin==0
    spikemaxwidth=7;
    updownmargin=0.4, 
    spikesign=0;
    [XX,FitX]=MakeDemoSpikeTrace;   
end

% ---------------------------------------------
 close all;
 Lx=length(XX);
 TT=(1:Lx)';
[step_props,~,~]=Get_StepsFromFit(TT,FitX);              

%% remove spikes; 
%run twice to clean ''comb sections''
     FitX_merged=FitX;
     for i=1:2
        [~,isolated,FitX_merged]=remove_spikes(step_props(:,1), step_props(:,5), spikemaxwidth, updownmargin, spikesign,FitX_merged);
        [step_props,~,~]=Get_StepsFromFit(TT,FitX_merged);
     end
    [StepFitProps_out,~,~]=Get_StepsFromFit(TT,FitX_merged);


if nargin <5
    plot(XX,'k-','LineWidth', 2); hold on;
    plot(FitX, 'c-','LineWidth', 1.5);
      plot(FitX_merged, 'r-', 'LineWidth', 2);
     legend('data','original fit','final fit');
     dum=1;
end

if nargin==3
     Axz=(1:length(XX))';
     data=[Axz XX FitX];
     outname1=strcat(PathName,'/',FileName(1:end-4),'_despiked.txt');
     dlmwrite(outname1,data);
     outname2=strcat(PathName,'/',FileName(1:end-4),'_despiked_props.txt');
     dlmwrite(outname2,StepFitProps_out);
end


    


 function  [spikes,isolated,FitX]=remove_spikes(index, steps, spikemaxwidth, updownmargin, spikesign,FitX);
%find all spikes; output is sets of indices:
%'spikes'  pairs of indices
%'isolated'

%indices and steps
OLi=index(1:end-3);     OLs=steps(1:end-3);
CLi=index(2:end-2);     CLs=steps(2:end-2);
CRi=index(3:end-1);     CRs=steps(3:end-1);
ORi=index(4:end);       ORs=steps(4:end);
%differences
dt_L=CLi-OLi;           ds_L=CLs-OLs;   %outer left, center left etc
dt_C=CRi-CLi;           ds_C=CLs-CRs;
dt_R=ORi-CRi;           ds_L=ORs-CRs; 

narrow=(dt_C<=spikemaxwidth);
opposed=(sign(CLs)+sign(CRs)==0);
nullingratio=abs((CRs+CLs))./(0.5*(abs(CRs)+abs(CLs)));
nulling=(nullingratio<updownmargin);
if spikesign~=0
    propersign=(sign(CLs)==spikesign);
else
    propersign=1+0*sign(CLs);
end

sel_L=find((narrow&opposed&nulling&propersign)==1)+1;  
 
spike_index_L=index(sel_L);
spike_index_R=index(sel_L+1);

spikes=sort(unique([spike_index_L spike_index_R]));
isolated=index(~ismember(index,spikes));

%remove the spikes
for spi=1:length(spike_index_L)    
    lft_i=spike_index_L(spi);
    rgt_i=spike_index_R(spi);
    padval=FitX(lft_i-1);
    FitX(lft_i:rgt_i)=padval;
end

%spike scanner:
%consider 4 points, 
    %-where the outer two are further away than the center;
    %where the two centers are:
    %shortly after one other
    %of approximate equal size
 %then, remove tese center two.
 %repeat this procedure to handle very dense blinking sections

dum=1;

function [XX,FitX]=MakeDemoSpikeTrace
    Lxi=2000;
    noize=10;                   %noise level
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
            %add spikes
            
            ups=round(linspace(105,Lxi-105,Nspikes));
            downs=ups+2;
            for spi=1:length(ups)
                peak=spike_ampli; %*randn(1,1);
                peaksign=2*round(rand(1,1))-1;
                XX(ups(spi):downs(spi))=XX(ups(spi):downs(spi))+peaksign*peak;
                FitX(ups(spi):downs(spi))=FitX(ups(spi):downs(spi))+peaksign*peak;
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