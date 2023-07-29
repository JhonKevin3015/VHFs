function [LXY1,LXY2,LYZ1,LYZ2,LXZ1,LXZ2,DRF,RF,LX,LY,LZ] = GRIDplanos(dirPC,dirdat)
%dirPC='/home/jhon/'
%dirdat=[dirPC 'MIT-JOB-TESIS/BAJARdatos/']  %  test
%dirdat='/Users/oceano/MIT-JOB-TESIS/run/' % para el grilladode 60X80
XC=rdmds([dirdat 'XC']);  % C para centrado
YC=rdmds([dirdat 'YC']);
DXC=rdmds([dirdat 'DXC']);
DYC=rdmds([dirdat 'DYC']);
XG=rdmds([dirdat 'XG']);   % geografico % en metros
YG=rdmds([dirdat 'YG']);
RC=rdmds([dirdat 'RC']);
%=========================
clear LX LY LZ
LX(:,1)=XG(:,1);
LY(:,1)=YG(1,:);
LZ(:,1)=RC(1,1,:);
%=========================
[LXY1 LXY2]=ndgrid(LX,LY);  % plano xy
[LXZ1 LXZ2]=ndgrid(LX,LZ); % plano xz
[LYZ1 LYZ2]=ndgrid(LY,LZ);  % plano yz
%=========================
DRF(:,1,1)=rdmds([dirdat 'DRF']);
RF(:,1,1)=rdmds([dirdat 'RF']);