%% ----------------
Initial Conditions
%% ----------------
clear all, close all
%---------------------
dirPC='/Volumes/JHON12T/';  
dirdatam='/Volumes/JHON12T/MIT-JOB-TESIS/model-dimentions/'
%===================================================
[LXY1,LXY2,LYZ1,LYZ2,LXZ1,LXZ2,DRF,RF,LX,LY,LZ] = GRIDplanos(dirPC,dirdatam);
[Y3 X3 Z3]=meshgrid(LY/1000,LX/1000,LZ);% con delta=3 se ve bien
%===================================================
nt=30;
dirdata=[../DATA/]
CAMPO='THETA';   
dat=rdmds([dirdata CAMPO],nt);
dat(:,1:2,:)=NaN;
dat(:,end-1:end,:)=NaN;
THETA=dat;
%=====
CAMPO='U';   
dat=rdmds([dirdata CAMPO],nt);
dat(:,1:2,:)=NaN;
dat(:,end-1:end,:)=NaN;
U=dat;
%======
[by,bx,bz,b] = fun_bonyance_and_gradient(dirdata,nt,LX,LY,LZ);
g=9.81;
rho0=1028.60;
RHO=-rho0*b/g+rho0;
ziseleter=30;
tAlpha=-2e-4   
dens=THETA*rho0*tAlpha+rho0;
%================================================