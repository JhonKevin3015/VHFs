addpath (genpath ([dirPC '../MATLABcodes/utils/matlab']))
addpath (genpath ([dirPC '../MATLABcodes/utils/matlab_funtions']))
%% =================================================================
%
%
%% ====================
clear all, close all
%======================
% All data are in DATA 
dirPC='../DATA/grid/';  
dirdatam='../DATA/grid/';
%===================================================
[LXY1,LXY2,LYZ1,LYZ2,LXZ1,LXZ2,DRF,RF,LX,LY,LZ] = GRIDplanos(dirPC,dirdatam);
[Y3 X3 Z3]=meshgrid(LY/1000,LX/1000,LZ);
%===================================================
dt=30
%========================
% activate this section to parallel option
%delete(gcp('nocreate'));
%parpool('local',4) 
%================================================
caso=1; % 2 or 3 or 4
dirdata=['../DATA/EXP';  num2str(caso) '/'];

nt=NT2(790) % chenges this     

%=================Partials==================
u=rdmds([dirdata 'U'],nt);
u(:,1:2,:)=NaN;
u(:,end-1:end,:)=NaN;
[uy ux uz]=gradient(u,LY,LX,LZ);
v=rdmds([dirdata 'V'],nt);
v(:,1:2,:)=NaN;
v(:,end-1:end,:)=NaN;
[vy vx vz]=gradient(v,LY,LX,LZ);
w=rdmds([dirdata 'W'],nt);
w(:,1:2,:)=NaN;
w(:,end-1:end,:)=NaN;       
[wy wx wz]=gradient(w,LY,LX,LZ);
%===============================================
[by,bx,bz,b] = fun_bonyance_and_gradient(dirdata,nt,LX,LY,LZ);
%===============================================
%%====================Physics parameters===========================
%=======Full Ertel-PV ===========================
f0=1e-4;   
q1=(wy-vz).*bx;
q2=(uz-wx).*by;
q3=(vx-uy+f0).*bz;
q=q1+q2+q3;
%==========%kinetic energy=====================
Ek=0.5*(u.^2+v.^2); 
%mask based in areas without motion =============
lin=reshape(Ek,[1000*300*60 1]);
lin(log2(lin)<-14)=NaN;
lin(isnan(lin)==0)=1;
MASK1=reshape(lin,[300 1000 60]);
MASK=Ek*0+1; % 
%=============================================
Div=ux+vy;  % divergencia
rv=vx-uy;   % vorticidad relativa
ss=vx+uy;   % 
sn=ux-vy;   %
S=sqrt(ss.^2+ sn.^2); % Deformation
OW=sn.^2+ss.^2-rv.^2; %Okubo-weiss
%============================================
ta=2*(-(ux.*bx + vx.*by)).*bx;
tb=2*(-(uy.*bx + vy.*by)).*by;
FF=ta+tb; % frontogenesis tendency
%===========================================
dt=30;
A=-2*(bx.^2 - by.^2).*sn*0.5;
B=-2*2*(bx.*by).*ss*0.5;
C=-2*(bx.^2 + by.^2).*Div*0.5;
Fsum=(A + B + C);tim(nn)=nt*dt/(60*60*24);
%============================================
%=================FF-front=====================
[Tprima] = fun_temp_anomaly(dirdata,dirPC,nt,caso);
cp=3850;ro=1025;
VHF=cp*ro*w.*Tprima;
%==============anomaly temp===============
%==============anomaly temp================================
nn=1
%================
pnz=25 ;%LZ(pnz) % from surface to depth corresponding to 100 m (LZ(25))
%================
y1=401 ; y2=600 ; %
%================
%OWm(nn,caso)=nanmean(OW(find(isnan(OW(:,:,1:pnz))==0)));
q=q.*MASK1; PVm(nn,caso)=nanmean(q(find(isnan(q(:,y1:y2,1:pnz))==0)));
FF=FF.*MASK1; FFm(nn,caso)=nanmean(FF(find(isnan(FF(:,y1:y2,1:pnz))==0)));
%Wm(nn,caso)=nanmean(W(find(isnan(W(:,:,1:pnz))==0)));
%Tprimam(nn,caso)=nanmean(Tprima(find(isnan(Tprima(:,:,1:pnz))==0)));
VHF=VHF.*MASK1; VHFm(nn,caso)=nanmean(VHF(find(isnan(VHF(:,y1:y2,1:pnz))==0)));
S=S.*MASK1; dati=S(:,y1:y2,1:pnz).^2;
    Sm(nn,caso)=sqrt(nanmean(dati(:)));
Div=Div.*MASK1; dati=Div(:,y1:y2,1:pnz).^2;
    Divm(nn,caso)=sqrt(nanmean(dati(:)));
Ekm(nn,caso)=nanmean(Ek(find(isnan(Ek(:,y1:y2,1:pnz))==0)));
%=================================================
%save ('TimeSerie.mat','PVm','FFm','VHFm','Sm','Divm','Ekm','tim')
% ==========================================================================


%   Figure time series  
%===============
load ('TimeSerie.mat')
%==========
close all
f1=figure('units','normalized','Position',[0.3 0.1 0.6 0.5]), hold on
%==========================================================================
%=========================================================
nc=1; nf=1 ; nnl=30;nnq=26
pos1 = [0.07+(0.28+0.055)*(nc-1) 0.57-(0.4+0.03)*(nf-1) 0.25 0.35];
    pl1=subplot('Position',pos1), hold on
    p2=plot(tim,Ekm(:,1),'--','color',rgb('black'),'LineWidth',3.5)
    p1=plot(tim,Ekm(:,4),'-','color','black','LineWidth',3.5)
    set(gca,'FontName','Iowan Old Style','FontSize',nnl);
    ylim([0.6 2.1]*10^(-3))
    ylabel('$\overline{\textbf{{KE}}}$ $\mathrm{(m^2/s^2)}$','FontSize',nnl,'Interpreter', 'latex')
    set(gca,'XTick',5:5:15,'XTickLabel',[])
    box on
    title('$\textbf{a)}$', 'Interpreter', 'latex','FontSize',32);
    
        hh=legend('Unforced','Forced');
        hh.FontSize=25
        hh.Position=[0.1 0.8 0.08 0.06]; 
        hh.Box='off'
%==========================================================================
%===========
nc=2; nf=1
pos1 = [0.07+(0.28+0.055)*(nc-1) 0.57-(0.4+0.03)*(nf-1) 0.25 0.35];
    pl1=subplot('Position',pos1), hold on
    plot(tim,Divm(:,1)*10^4 ,'--','color',rgb('black'),'LineWidth',3.5)
    plot(tim,Divm(:,4)*10^4,'-','color','black','LineWidth',3.5)
    set(gca,'FontName','Iowan Old Style','FontSize',nnl); 
    ylim([0 10]*10^-1)
    ylabel('$\textbf{RMS}$ ($\mathrm{\delta /f}$)','FontSize',nnl,'Interpreter', 'latex')
    set(gca,'XTick',5:5:15,'XTickLabel',[])
    box on
    title('$\textbf{b)}$', 'Interpreter', 'latex','FontSize',32);

%==========================================================================
%===========
nc=3; nf=1
pos1 = [0.07+(0.28+0.055)*(nc-1) 0.57-(0.4+0.03)*(nf-1) 0.25 0.35];
    pl1=subplot('Position',pos1), hold on
    plot(tim,Sm(:,1)*10^4,'--','color',rgb('black'),'LineWidth',3.5)
    plot(tim,Sm(:,4)*10^4,'-','color','black','LineWidth',3.5)
    set(gca,'FontName','Iowan Old Style','FontSize',nnl);
    ylim([0 10]*10^-1)
    ylabel('$\textbf{RMS}$ ($\mathrm{S/f}$)','FontSize',nnl,'Interpreter', 'latex')
    set(gca,'XTick',5:5:15,'XTickLabel',[])
    title('$\textbf{c)}$', 'Interpreter', 'latex','FontSize',32);
    box on  
%==========================================================================
%===========
%==========================================================================
nc=1; nf=2
pos1 = [0.07+(0.28+0.055)*(nc-1) 0.55-(0.4+0.03)*(nf-1) 0.25 0.35];
    pl1=subplot('Position',pos1), hold on
    plot(tim,PVm(:,1),'--','color',rgb('black'),'LineWidth',3.5)
    plot(tim,PVm(:,4),'-','color','black','LineWidth',3.5)
    set(gca,'FontName','Iowan Old Style','FontSize',nnl);
    ylim([-2 3.5]*10^-10)
    ylabel('$\overline{\textbf{q}}$ $\mathrm{(s^{-3})}$','FontSize',nnl,'Interpreter', 'latex')
    set(gca,'XTick',5:5:15,'XTickLabel',[5 10 15])
    title('$\textbf{d)}$', 'Interpreter', 'latex','FontSize',28);
    box on        
%==========================================================================
nc=2; nf=2
pos1 = [0.07+(0.28+0.055)*(nc-1) 0.55-(0.4+0.03)*(nf-1) 0.25 0.35];
    pl1=subplot('Position',pos1), hold on
    plot(tim,FFm(:,1),'--','color',rgb('black'),'LineWidth',3.5)
    plot(tim,FFm(:,4),'-','color','black','LineWidth',3.5)
    set(gca,'FontName','Iowan Old Style','FontSize',nnl);
    ylim([-0.02 1.0]*10^-18)
    ylabel('$\overline{\textbf{Fs}}$ $\mathrm{(s^{-5})}$','FontSize',nnl,'Interpreter', 'latex')
    box on
    xlabel('$\textbf{{Time}}$ $\mathrm{(days)}$','FontSize',nnl,'Interpreter', 'latex')
    set(gca,'XTick',5:5:15,'XTickLabel',[5 10 15])
    title('$\textbf{e)}$', 'Interpreter', 'latex','FontSize',32);
%===============================================================
%===============================================================
nc=3; nf=2 ; cp=3850;ro=1025
pos1 = [0.07+(0.28+0.055)*(nc-1) 0.55-(0.4+0.03)*(nf-1) 0.25 0.35];
    pl1=subplot('Position',pos1), hold on
    plot(tim,VHFm(:,1),'--','color',rgb('black'),'LineWidth',3.5)
    plot(tim,VHFm(:,4),'-','color','black','LineWidth',3.5)
    set(gca,'FontName','Iowan Old Style','FontSize',nnl);
    ylim([-0.5 35])
    ylabel('$\overline{\textbf{VHF}}$ $\mathrm{(W^2/m^2)}$','FontSize',nnl,'Interpreter', 'latex')
    set(gca,'XTick',5:5:15,'XTickLabel',[5 10 15])
    title('$\textbf{f)}$', 'Interpreter', 'latex','FontSize',32);
    box on 
%==========================================================================
set(gcf,'Color','white','Renderer','zbuffer');
set(gca,'Color','white','XColor','black', ...
    'YColor','black');
%===================================================================