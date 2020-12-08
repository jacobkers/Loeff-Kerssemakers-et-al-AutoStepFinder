function ProbeStep_Extract(tol_idx,tol_step);
% description: run after using ProbeStepInject and Autostepfinder
% ProbeStepExtract(margin_idx)
% 
% parameters: margin_idx is the tolerance, in time points, that a step can
% be off-place to be accepted. default is 0 (no offset)

%output: percentage of correct placements; the error is based on the
%variation per trace
% 
%example: ProbeStepExtract(5,7) returns the percentage of injected steps that were
%found correctly within 5 time points and 7 signal units of the correct value 

%Reference: Auxiliairy tool of AutoStepfinder: 
%A fast and automated step detection method for single-molecule analysis.
%Luuk Loeff*, Jacob Kerssemakers*, Chirlmin Joo & Cees Dekker.
% * Equal contribution
% Last update by Jacob Kerssemakers, October 2020
%--------------------------------------------------------------------------
if nargin==0,
    tol_idx=10;
    tol_step=100;
end

codepth=pwd;  
mpath=pwd; %you may replace this by another path
infopath=strcat(mpath, '\ProbeStepInject\'); 
fitpath=strcat(infopath, 'traces\StepFit_Result\'); 

bins=25;
cd(fitpath); FilNames=dir('*_fits.txt'); cd(codepth); croptextno=9;

%collect the list of original steps 
StepInjectPropsName='probe_step_info.txt';   
oriprops=double(dlmread(strcat(infopath,StepInjectPropsName)));
%[repeat index stepsize];
Nrep=length(FilNames);
for rp=1:Nrep
    %get numeric repeat number
    repeatfilename=FilNames(rp).name;
    repno_i = strfind(repeatfilename,'repeat')+6;
    repno=str2num(repeatfilename(repno_i:repno_i+3));
    sel=find(oriprops(:,1)==rp);
    ori_id=oriprops(sel,2);  %original indices
    ori_st=oriprops(sel,3);  %original steps
    
    fitpropsname=strcat(FilNames(rp).name(1:end-croptextno),'_properties.txt');
    fitprops=double(dlmread(strcat(fitpath,fitpropsname),',',1,0));
    fit_id=fitprops(:,1);   %fit indices
    fit_st=fitprops(:,5);   %fit steps
    %check indices 
    index_tolerance=tol_idx/max([ori_id ; fit_id]);
    [somethingnearby,which_ones_in_fit]=ismembertol(ori_id,fit_id, index_tolerance);
    found_ori_idx=(ori_id(somethingnearby==1));
    %check stepsizes of these   
    found_ori_st=ori_st(somethingnearby==1);                    %the steps that were found
    which_idx=sort(which_ones_in_fit((which_ones_in_fit>0)));  %corresponding fit locations
    candidate_st=fit_st(which_idx);                             %corresponding step sizes
    stepsizedif=abs(found_ori_st-candidate_st);
    
    found_nearby_at_similarsize=find(stepsizedif<tol_step);
    
    
    gooddetection=length(found_nearby_at_similarsize);
    
    perc(rp)=100*sum(gooddetection)/length(ori_id);
    dum=1;
end
probe_perc=mean(perc);
st_perc=std(perc);
err_perc=st_perc/(Nrep^0.5);
disp(['detection percentage mean+SEM=',num2str(probe_perc, '%02.1f'),'%'...
     ' +/- ' num2str(err_perc, '%02.1f')]);
    
    dum=1;