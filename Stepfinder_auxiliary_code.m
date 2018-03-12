function Stepfinder_auxiliary_code
%% Various useful functions to be used with the main Stepfinder code for Matlab users






%% This section (xxx to ~xxxx) contains the 'Core' function of the stepfinder; 
%it can be cut and autorun independently (on an auto-simulated curve) for demo purposes
function [FitX,stepsX,S_fin]=StepfinderCoreCopy(X,initval)
%This function performs a step-fitting algorithm on data in a quick fashion.
%Background theory can be found in  Supplement 3 of the paper: "Assembly dynamics of microtubules at molecular resolution" 
%J.W. Kerssemakers, E. L. Munteanu, L. Laan, T.L. Noetzel, M.E. Janson, M. Dogterom 
%Nature  442(7103) 709-12 (2006). 
%Algorhitm code developed by Jacob Kerssemakers
%GUI shell developed by Luuk Loeff
%output: list of stepsizes: [index time  levelbefore levelafter step dwelltimeafter steperror]
if nargin<2
    close all;
    addpath('D:\jkerssemakers\My Documents\BN CD Recent\BN_CD16 LuukLoeff\HandyTools')
    X=BottomUpTraceBuilder;
    initval.stepnumber=floor(length(X)); 
    initval.overshoot=1;
    initval.DemoMode=1;
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
if initval.DemoMode
    close all;
    subplot(2,1,1); plot(X,'y'),hold;plot(FitX,'k','LineWidth',2);
    title('Data and Fit');xlabel('time');ylabel('position,a.u.');
    subplot(2,1,2); semilogx(S_fin,'k-','LineWidth',2);
    title('S-curve');xlabel('Stepnumber');ylabel('S-value, a.u.');
end


function [bestshot,S_fin]=Eval_Scurve(S_raw)
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
        
        plotmoment=((c<10)|((c>50)&(c<1000)&(mod(c,10)==0))|...
                    (c>1000)&(c<2000)&(mod(c,100)==0)|...
                    (c>2000)&(mod(c,500)==0));
                
        if initval.DemoMode && plotmoment
            %close all;
            subplot(2,1,1); 
                plot(X,'y'),hold on;
                plot(FitX,'k','LineWidth',1);
                hold off;
            title('Data and Fit');xlabel('time');ylabel('position,a.u.');
            subplot(2,1,2); 
                semilogx(S,'-k','LineWidth',2);
                title('S-curve');xlabel('Stepnumber');ylabel('S-value, a.u.');
                hold off;
            pause(0.01);
        end
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
    
    
function xx=BottomUpTraceBuilder
%JWJK_B:-------------------------------------------------------------------
%Self-scaling step trace
%
%Summary: this function builds self-scaling step curves to illustrate
%properties of 'S-curve' based stepfinding
%
%Approach: A self-scaling step trace is built as follows: We start with two points
%separated by unit step size, i.e. a unit step between two plateasu of just one point. 
%We concatenating this section Np times with copies of itself. 
%This creates a plateau of Np*2 points. Next,this procedure is repeated on the plateau , 
%now with step size two. In total, this concatenation step is repeated Nc times, 
%increasing the step size from 1 to Nc. Finally, this creates a 'fractal' 
%curve of 2*(Np*2)^Nc points, resembling a square wave lines with smaller 
%square waves which in turn are %lined with even smaller square waves, and so on. 
%Such a curve contains steps at Nc-1 length scales, 
%resulting in as any peaks in the S-curve. The relative abundance and size
%of per 'step scale' detrmines which peak will dominate, 
%since a peak height scales linearly with the abundance of steps and quadratically 
%with the size of associated steps. For  example, with a step size increase of 1 
%and a number increase of factor 2*Nc, peaks will roughly be equal. 
%Note that the rising 'noise tail' in the S-curve, as observed for experimental data, 
%here amounts to the final fitting of unit step sizes and can thus be
%interpred as the 'Nc_th S-peak'.
 
%Jacob Kers 2017
%
%Input: repaat factors
%
%Output: trace
%
%References: Jacob Kers 2017
%
%:JWJK_B-------------------------------------------------------------------

    Nc=4; %Concatenation doubling repeats (amount to number of S_peaks+1);
    Np=5; %plateau repeats (pushes 'noise tail' to higher Nsteps)
    StepN=1; %increase of step size per 
    xx=[-0.5 0.5] ; %the start unit
    for nn=1:Nc
        stepN=StepN*(2^nn); %steps are growing to keep S-function peaks flat
        leftpart=xx-stepN/2;
        rightpart=xx+stepN/2;
        buf=[leftpart rightpart];
        xx=[[repmat(buf,1,Np)]];
    end   
    
    