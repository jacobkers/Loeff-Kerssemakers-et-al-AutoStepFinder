function ProbeStep_Inject(stepspacing,stepshifts, stepsize,datain)
%JWJK_B:-------------------------------------------------------------------
% data_out=ProbeStepInject(spacing,shift, size, data_in)  
% 
%use:
%1) run first this tool; traces are saved in a local directory '\ProbeStepInject\traces' 
%2  batch-run Autostepfinder on this directory
%3) run ProbeStepExtract (see settings there)
%
% description: 
% ProbeStepInject(stepspacing,stepshift, stepsize, data_in)  uses data
% ProbeStepInject(stepspacing,stepshift, stepsize)  user loads data from pwd
% ProbeStepInject() makes demo data
% 
% parameters:
% spacing: the plateau in pts between up-and-down steps.
% [note on density]
% shift: index where the step series starts. 
% To boost statistics, shift can be a vector; then multiple saves will be done, each with a shift specified by this vector. 
% size: up-and-down step size in units of the input data
% data_in: single-column data
% 
%example: ProbeStepInject(100,1:10:100, 10) user-loads a trace injects a blockwave with
%plateaulength 100, amplitude 10, and repeats the action ~10 times with
%increasing placement shift.

%Output: saves one or multiple traces to be analyzed by AS. As an example,
%the first trace is shown.
%



%Extras: includes autorun demo option of a  noisy step trace

%Reference: Auxiliairy tool of AutoStepfinder: 
%A fast and automated step detection method for single-molecule analysis.
%Luuk Loeff*, Jacob Kerssemakers*, Chirlmin Joo & Cees Dekker.
% * Equal contribution
% Last update by Jacob Kerssemakers, October 2020
%--------------------------------------------------------------------------

%:JWJK_B-------------------------------------------------------------------
shoit=0;
if nargin==0
    close all;
    shoit=1;
    stepsize=10;
    stepspacing=500; 
    stepshifts=0:10:90;
    datain=MakeDemoTrace;
    Label=strcat('demo');
end

if (nargin==3)  %Load
    shoit=1;
    [FileName,PathName] = uigetfile('.txt', 'Open the file of interest');
    datain=double(dlmread(strcat(PathName,'/',FileName)));
    [pts,cols]=size(datain);
    if cols>1, trace=datain(:,2); else trace=datain;end   
end

outdir=strcat(pwd,'\ProbeStepInject\');
if isdir(outdir), rmdir(outdir,'s'),end;
mkdir(outdir);
mkdir(outdir,'traces')
allinfo=[];

for rp=1:length(stepshifts)
    disp(['saving trace:',num2str(rp),'of',num2str(length(stepshifts))]);
stepshift=stepshifts(rp);
stepsize=10;
L_total=length(datain);
L_quarterblock=stepspacing/2;
reps=ceil(L_total/(4*L_quarterblock));
      
blocktemplate=[zeros(1,L_quarterblock) ones(1,2*L_quarterblock)...
               zeros(1,L_quarterblock)];
blockcurve=[];
for ii=1:reps
    block=stepsize*blocktemplate;
    blockcurve=[blockcurve block];
end


%shifting, cropping
if stepshift>0
    addit=zeros(1,stepshift)+blockcurve(1);
    blockcurve=[addit blockcurve];
end
blockcurve=blockcurve(1:L_total);

LB=length(blockcurve);
if LB<L_total&&LB>0  %padding
    blockcurve=[blockcurve zeros(1,L_total-LB)+blockcurve(end)];
end



%% measure the location of injected steps
if LB>0
    steplocs=diff(blockcurve);
    indexes=(find(steplocs~=0))';
    steps=(steplocs(indexes))';
    reps=0*steps+rp;
    dataout=(datain'+blockcurve)';
else
    indexes=NaN;
    steps=NaN;
    reps=0*steps+rp;
    blockcurve=0*datain';
    dataout=(datain+blockcurve');
end

allinfo=[allinfo; [reps indexes steps]];


if shoit&&rp==1
    close all;
    plot(datain,'k-'); hold on;
    plot(blockcurve-0.5*stepsize,'b-'); hold on
    plot(datain+blockcurve'+4*stepsize,'r-');  
    legend('original','probing steps','summed' );
    title('spiking of real data with flat-range, equidistant steps');
    xlabel('time, a.u.');
    ylabel('extension, nm');
    pause(0.5);
end

%% save the curve


outname=strcat('test_repeat',num2str(rp, '%04.0f'),'.txt');
dlmwrite([outdir,'traces\',outname],dataout);

end
stepinfo=[indexes steps];
outname2='probe_step_info.txt';
dlmwrite([outdir,outname2],allinfo);

dum=1;

function XX=MakeDemoTrace
    Lxi=20000;
    noize=5;                   %noise level
    XX=[]; FitX=[];
    stepsizes=[-10 10 30 50 5];
    Nspikes=30; 
    spike_ampli=100;
    isolated=[1 [10:500:Lxi]];  %just a step train
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
            end 

