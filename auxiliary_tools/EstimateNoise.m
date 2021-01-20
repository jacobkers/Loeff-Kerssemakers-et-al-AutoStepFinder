function est_noise=EstimateNoise(trace)
%% Get an estimation of noise in a trace
%Aim: get broad estimate of noise levels in a stepped signal.
%Approach: pairwise differences are taken between time points to yield a
%global noise estimate. The distance between time points is varied. A plot
%of noise level as function of time delay (in points) is shown to the user; 
%user can indicate a representative noise
%level

%Input: 
%1) EstimateNoise: no input, runs a demo trace
%2) EstimateNoise(1): load data as .tx column file. For more than one column, 
    %we assume %column 1 is time and only column 2 is processed.
%3) EstimateNoise(trace)
    %trace is a single or multiple data column as in (2). 

%Output: user-estimate of noise level

%Extras: includes autorun demo option of a low-pass fltered noisy step trace

%Reference: Auxiliairy tool of AutoStepfinder: 
%A fast and automated step detection method for single-molecule analysis.
%Luuk Loeff*, Jacob Kerssemakers*, Chirlmin Joo & Cees Dekker.
% * Equal contribution
% Last update by Jacob Kerssemakers, October 2020
%--------------------------------------------------------------------------

tracelimit=10000;   %crop long traces
max_range=100;      %maximum disance between time points to consider


%% Get or build trace
if (nargin==1)  %Load
    [FileName,PathName] = uigetfile('.txt', 'Open the file of interest');
    trace0=double(dlmread(strcat(PathName,'/',FileName)));
    [pts,cols]=size(trace0);
    if cols>1, trace=trace0(:,2); else trace=trace0;end   
end
%Build trace
if nargin < 1
    close all; 
    LT=2000;
    N0=20;
    dstep=(20+10*rand(20,1));
    istep=cumsum(round(100+00*rand(N0,1)));
    trace=10*randn(1,LT);
    for jj=1:N0
        trace(istep(jj):end)=trace(istep(jj):end)+dstep(jj);
    end
   %trace=smooth(trace',6)';
   trace=trace';
end

%% Get pairwise difference as function of distance ------------------------

%crop
if length(trace)>tracelimit, trace=trace(1:tracelimit),end
Cp=length(trace);

%build difference maps
map=repmat(trace',max_range,1);
shiftmap=NaN*zeros(max_range,Cp);
for ii=1:max_range
    shiftmap(ii,1:Cp-ii)=trace(ii+1:end);
end
difmap_sq=(map-shiftmap).^2;

%get differences
crop_hi=0.01;
noisecurve=zeros(max_range,1);
for rg=1:max_range
    sqdif_sort=sort(difmap_sq(rg,:));
    sqdif_lo=sqdif_sort(1:round((1-crop_hi)*Cp));
    noisecurve(rg)=(nanmean(sqdif_lo)).^0.5/(2^0.5);
end

noise_nn=0*noisecurve+noisecurve(1);

%% Plot and user interaction menu
close all;
figure('Name','Noise_Estimator','NumberTitle','off','units', 'normalized', 'position', [0.745 0.1 0.25 0.4]);
  %[0.745 0.32 0.25 0.6]);
    plot(noisecurve,'LineWidth',2, 'Color',[0,0.2,1]);  hold on
    plot(noise_nn,'LineWidth',1, 'Color','k');  hold on       
    xlabel('time difference, pts')
    ylabel('noise, measurmeent units');
    ylim([0 max(noisecurve)]);
    legend('noise vs delay','nearest-neighbour noise', 'Location', 'SouthOutSide');
    title('Click for user estimate')
    [~,est_noise,~]=ginput(1);
    noisecurve_est=0*noisecurve+est_noise;
    tx=text(max_range/10,est_noise+1,num2str(est_noise,'% 3.2f'));
    tx(1).Color = 'red';
    tx(1).FontSize = 12;
    plot(noisecurve_est,'LineWidth',1, 'Color','r');  hold on
    title('Including user estimate')
    legend('noise vs delay','nearest-neighbour noise', 'user-estimate','Location', 'SouthOutSide');
    
