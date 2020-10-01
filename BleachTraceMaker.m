function MakeBleachTrace
%Make a simple bleach-like trace

pts=1000; tau=pts/10; N0=20; noise=0.5;

tt=linspace(1,pts,1000);
yy=round(N0*exp(-tt/tau))+noise*randn(1,pts);

close all;
plot(tt,yy);
data=[tt', yy'];

SaveName=strcat('testdata_simulated_BleachTrace',num2str(N0),'steps.txt');
dlmwrite(strcat(pwd,'/',SaveName), data);
