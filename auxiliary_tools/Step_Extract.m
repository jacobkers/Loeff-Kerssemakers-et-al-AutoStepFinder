function Step_Extract
%JWJK_B:-------------------------------------------------------------------
%Step Extraction
%
%Summary: This tool processes the saved results of a step finding analysis.
%
%Input: step finder output files
%
%Output: plots with histograms
%
%References: 
%[1]  Real-time detection of condensin-driven DNA compaction reveals a multistep binding mechanism
% Authors, Jorine Mirjam Eeftens, Shveta Bisht, Jacob Kerssemakers, Christian Haering, Cees Dekker
% Publication date 2017/1/1, Journal bioRxiv
% [2] Assembly dynamics of microtubules at molecular resolution.
% Nature. 2006 Aug 10;442(7103):709-12. Epub 2006 Jun 25.
% Kerssemakers JW1, Munteanu EL, Laan L, Noetzel TL, Janson ME, Dogterom M.
% [3] Loeff, Kerssemakers et al 2017
%:JWJK_B-------------------------------------------------------------------

mpath=pwd;
mpath='C:\Users\jkerssemakers\CD_Data_in\2019_Richard\20191029_pilots\Trimmed traces 1 nM Condensin 50 uM ATP 0.3 pN_adapt\';
codepath=pwd;  
fitpath=strcat(mpath, 'StepFit_Result\');  
bins=40;
batchload=1;
 
croptextno=9;
if batchload
    cd(fitpath)
    FilNames=dir('*_fits.txt'); cd(codepath); 
else
    fitpath=pwd;
    cd(fitpath)
    [FileName,PathName] = uigetfile('*_fits.txt','Select the "_fits" file');
    FilNames.name=FileName;
    FilNames.folder=PathName;
    fitpath=strcat(PathName);
end
cd(codepath);

AllFiles=length(FilNames);
%FilNames(1).name

AllIndex=[];
AllTime=[];
AllLevelBefore=[];
AllLevelAfter=[];
AllLevelMid=[];
AllStepSizes=[];
AllDwellTimeBefore=[];
AllDwellTimeAfter=[];
AllMeasuredError=[];        
AllNoise=zeros(AllFiles,1);
AllS1max=zeros(AllFiles,1);
AllS2max=zeros(AllFiles,1);

close all

for ii=1:AllFiles
AllFiles-ii+1
%read in all the fit information-------------------------------------------

StepFitsName=strcat(FilNames(ii).name(1:end-croptextno),'_fits.txt');
StepFits=double(dlmread(strcat(fitpath,StepFitsName),',',1,0));       
StepPropsName=strcat(FilNames(ii).name(1:end-croptextno),'_properties.txt');
StepFitProps=double(dlmread(strcat(fitpath,StepPropsName),',',1,0));
SCurveName=strcat(FilNames(ii).name(1:end-croptextno),'_s_curve.txt');
SCurveProps=double(dlmread(strcat(fitpath,SCurveName),',',1,0));


AllS1max(ii)=SCurveProps(1,2);
AllS2max(ii)=SCurveProps(2,2);


%Traces: Columns Time ,Data, FinalFit
TT=StepFits(:,1);    %time axis
XX=StepFits(:,2);    %raw data (on which step fit was done)
FitXX=StepFits(:,3); %step fit
Residu=XX-FitXX;
AllNoise(ii)=std(Residu);


%Step Properties: 
%IndexStep,TimeStep,LevelBefore,LevelAfter,StepSize,
%DwellTimeStepBefore,DwellTimeStepAfter,StepError

Index=StepFitProps(:,1);     
Time=StepFitProps(:,2);    
LevelBefore=StepFitProps(:,3);     
LevelAfter=StepFitProps(:,4); 
LevelMid=(LevelBefore+LevelAfter)/2;
StepSizes=StepFitProps(:,5);    
DwellTimeBefore=StepFitProps(:,6);     
DwellTimeAfter=StepFitProps(:,7);
MeasuredError=StepFitProps(:,8);

%collect them
AllIndex=[AllIndex; Index];
AllTime=[AllTime; Time];
AllLevelBefore=[AllLevelBefore; LevelBefore];
AllLevelAfter=[AllLevelAfter; LevelAfter ];
AllLevelMid=[AllLevelMid; LevelMid];
AllStepSizes=[AllStepSizes; StepSizes];
AllDwellTimeBefore=[AllDwellTimeBefore; DwellTimeBefore];
AllDwellTimeAfter=[AllDwellTimeAfter; DwellTimeAfter];
AllMeasuredError=[AllMeasuredError; MeasuredError];    


    
if 0 
close all;
disp(FilNames(ii).name);
figure(1); 
plot(Residu);
hold on;

title(FilNames(ii).name);
xlabel('time');
%ylabel('position,a.u.');
[~]=ginput(1);
pause(0.5); 
end        
end

%% histogram all steps--------------------------------------------- 
    minstepsize=nanmin(AllStepSizes);
    maxstepsize=nanmax(AllStepSizes);
    symXplotrange=1.1*max([abs(minstepsize) abs(maxstepsize)]);
    hx=linspace(minstepsize,maxstepsize,bins);
    HistX=(hist(AllStepSizes,hx)'); 
    
%% histogram all dwells--------------------------------------------- 
    mindwelltime=nanmin(AllDwellTimeAfter);
    maxdwelltime=nanmax(AllDwellTimeAfter);
    Tplotrange=1.1*max([abs(mindwelltime) abs(maxdwelltime)]);
    ht=linspace(mindwelltime,maxdwelltime,bins);
    HistT=(hist(AllDwellTimeAfter,ht)'); 


%% plotting panels
subplot(2,2,1);
plot(TT,XX,...                                                  %Plot Data
        'LineWidth',2,....                                                  %Linewidth
        'Color',[0,0.2,1]);                                                 %Color line RBG
        hold on
        plot(TT,FitXX,...                                              %Plot Fit
        'LineWidth',1,....                                                  %Linewidth
        'Color',[1,0.7,0]);                                                 %Color line RBG
        title('Last fitted trace');
        MaxX=Time(end);                                             %Determine length X axis
        MaxY=max(XX)*1.2;                                         %Determine length Y axis
        MinY=min(XX);                                         %Determine length Y axis
        xlim([0 MaxX]);                                             %Set X axis
        ylim([MinY MaxY]);                                  %Set Y axis
        xlabel('Time (s)','FontSize',12);                                   %Label X axis
        ylabel('Position (A.U.)','FontSize',12);                            %Label Y axis

 subplot(2,2,2);
plot(SCurveProps(:,1),SCurveProps(:,2:3));

 subplot(2,2,3);   
    bar(hx,HistX,'y');        
    title('Step Histogram');
    xlabel('Step size, a.u.'); ylabel('Counts');
    xlim([-symXplotrange symXplotrange]);        

    subplot(2,2,4);   
    bar(ht,HistT,'b');        
    title('Dwell Histogram');
    xlabel('Dwell Time, a.u.'); ylabel('Counts');
    xlim([0 Tplotrange]);      

disp(strcat(num2str(length(AllIndex)),'steps'));