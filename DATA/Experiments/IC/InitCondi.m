%===================INITIAL CONDITIONS FOR EXPERIMENTS ====================
clear all, close all
ny=1000;%360; long-channel 
nx=300;%360; cross-channel
nz=60; % depth 
deltaxy=200%200 % metros
perturba=1 % 1 Add white noise into the density field
%==========================================================================
[dens,theta,uinit,etainit,zc,dh,Us,Ri,Ls,Lsmin,ts] = frontInitCondi(ny,nx,nz,deltaxy,perturba);
%==========================================================================
tem3dN=theta;  % temp
uvelN=uinit;      % zonal velocity
ETA2dN=etainit; % surface lavel 
%================================
SA3dN=tem3dN*0; 
vvelN=uvelN*0;
x=[1:nx]; y=[1:ny]; 
z=-zc;
delR = dh;
%----------------------------------------------
dirsave='YourPath'
Savebinary(dirsave,tem3dN,ETA2dN,SA3dN+34,uvelN,vvelN,x,y,z,delR)
%==========================================================================
%==============Qnet: Forcing===========================
            HORA=[0:0.5:hora(end)]
            QnetIDE=(-cos(2*pi*HORA/24)+1)*150
            QnetIDE(1,1:2*24*2)=0;
            figure (); plot(HORA/24, QnetIDE), xlabel('time (days)'), ylabel('Qnet'), ylim([-5 305])
            %===============================================================
            for nt=1:length(HORA)
                    QnetSave(1:nx,1:ny,nt)=QnetIDE(nt);
            end
%=======================================================
fid=fopen([dirsave 'QnetIDE.forcing'],'w','b');
fwrite(fid,QnetSave,'float32');fclose(fid);
%===================================================================
