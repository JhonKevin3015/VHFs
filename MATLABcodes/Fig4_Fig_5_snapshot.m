addpath (genpath ([dirPC '../MATLABcodes/utils/matlab']))
addpath (genpath ([dirPC '../MATLABcodes/utils/matlab_funtions']))
%% ----------------
%Initial Conditions
%Uses the Outputs of the simulation to time=30 s. It is close to time=0. 
%% ----------------
clear all, close all
%-------------------
% All data are in DATA  
dirPC='../DATA/grid/';  
dirdatam='../DATA/grid/';
%===================================================
[LXY1,LXY2,LYZ1,LYZ2,LXZ1,LXZ2,DRF,RF,LX,LY,LZ] = GRIDplanos(dirPC,dirdatam);
[Y3 X3 Z3]=meshgrid(LY/1000,LX/1000,LZ);% 
%===================================================
%================================================
caso=1; %Unforced case with KPP (Fig4)
%=========================
dirdata=['../DATA/EXP';  num2str(caso) '/'];
NT=30:60:(960*60)-30;%NT=30:60:(1440*60)-30;
%=============================================================
close all
f1=figure('units','normalized','Position',[0.3 0.1 0.65 0.65]), hold on
%=====================================================================        
con=0
for nt=NT(700)%
%==============vertical velocity===============
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
%==============anomaly temp===============
[Tprima] = fun_temp_anomaly(dirdata,dirPC,nt,caso);
%=================Okubo-weiss==================
ta=2*(-(ux.*bx + vx.*by)).*bx;
tb=2*(-(uy.*bx + vy.*by)).*by;
FF=ta+tb; % frontogenesis tendency
f0=1e-4;  
%=======PV-ERtel=================================
q1=(wy-vz).*bx;
q2=(uz-wx).*by;
q3=(vx-uy+f0).*bz;
q=q1+q2+q3;
%================================================
Div=ux+vy;  % Divergence
rv=vx-uy;   % Relative Vorticity
ss=vx+uy;   % 
sn=ux-vy;   %
S=sqrt(ss.^2+ sn.^2); % Deformation
OW=sn.^2+ss.^2-rv.^2; %Okubo-weiss
cp=3850;ro=1025;
VHF=cp*ro*w.*Tprima; % vertical heat flux
%================================================
aa=34;bb=26;
%f1=figure('units','normalized','Position',[0.3 0.1 0.65 0.65]), hold on
nz=14; % 50 meters deep.
nc=3; nf=1
pos1 = [0.07+(0.28+0.02)*(nc-1) 0.55-(0.4+0.03)*(nf-1) 0.28 0.4];
    pl1=subplot('Position',pos1)
        hold on
        %============PV========
        pcolor(LXY1/1000, LXY2/1000-100,q(:,:,nz)); shading interp;    
        %[hh cc]=contour(LXY1/1000, LXY2/1000,PV(:,:,nz),[-.1,-.1]*10^(-9),'ShowText','off','color',rgb('green'));
        %cc.LineWidth=1.5, cc.LineStyle='-'
        [hh cc]=contour(LXY1/1000, LXY2/1000-100,q(:,:,nz),[0,0],'ShowText','off','color',[0.6 0.6 0.6]);
        cc.LineWidth=0.8, cc.LineStyle='-' ;  
        caxis([-2 2]*10^(-9));
        colormap(BWR2(360));
        %hh=colorbar
        %cmocean('balance','pivot',0)
    %    title('Ertel-PV ($1/s^3$)', 'Interpreter','latex','FontSize',bb), %colorbar ;
        set(gca,'FontName','Iowan Old Style','FontSize',aa);
        axis equal,xlim([0 60]), ylim([-20 20])
        set(gca, 'box','on','XTickLabel',[ 9 18 27 ],'XTick',[],'YTickLabel',[9 18 27 ]-18,'YTick',[]-18)
        h3=text(1,21,{'$\textbf{c)}$'},'color','black',...
    'FontSize',35.0,'HorizontalAlignment','left','VerticalAlignment','bottom','Interpreter','latex','Rotation',0);
%=====================
    h1=plot([0 60 60 0 0] ,[-20 -20 20 20 -20],'-','color','black', 'LineWidth',1.5)       
%=====================
    hh=colorbar(pl1,'Ticks',[-2:2:2]*10^(-9) ,'TickLabels',...
    {'$-2$' '$0$' '$2$'},'Location','northoutside',...
    'Position',[ 0.07+(0.28+0.02)*(nc-1)+0.15 0.55-(0.4+0.06)*(nf-1)+0.30  0.1  0.01],'TickLabelInterpreter','latex','FontSize',32.0); 
    set(get(hh,'title'),'string',[{'$\mathrm{q (s^{-3}) \times 10^{-9}}$'}],'FontSize',32.0,'Interpreter','latex','color',rgb('black')); 
    hh.Box= 'on'
    hh.Color= rgb('black')  
    colormap(BWR2(360));
%==========================================================================    
nc=1; nf=1
pos1 = [0.07+(0.28+0.02)*(nc-1) 0.55-(0.4+0.03)*(nf-1) 0.28 0.4];
    pl1=subplot('Position',pos1), hold on
    pcolor(LXY1/1000, LXY2/1000-100,OW(:,:,nz)), shading interp
    %[hh cc]=contour(LXY1/1000, LXY2/1000,PV(:,:,nz),[-.1,-.1]*10^(-9),'ShowText','off','color',rgb('green'));
    %cc.LineWidth=1.5, cc.LineStyle='-'
    caxis([-3 3]*10^(-8))
    colormap(BWR2(360))
    %cmocean('balance','pivot',0)
    %title('Okubo-Weiss ($1/s^2$)', 'Interpreter', 'latex','FontSize',bb)% colorbar 
    set(gca,'FontName','Iowan Old Style','FontSize',aa);
        axis equal,xlim([0 60]), ylim([-20 20])
    set(gca, 'box','on','XTickLabel',[ 9 18 27 ],'XTick',[],'YTickLabel',[-15 0 15],'YTick',[-15 0 15])
    ylabel('Y-axis (km)', 'Interpreter', 'latex')%,ylabel('Profundidad (m)', 'Interpreter', 'latex');
    dt=30
            h3=text(1,21,{'$\textbf{a)}$'},'color','black',...
    'FontSize',35.0,'HorizontalAlignment','left','VerticalAlignment','bottom','Interpreter','latex','Rotation',0);
%    text(18,-5,['Time: ' num2str(round(nt*dt/(60*60),2)) ' hours (' num2str(round(nt*dt/(60*60*24),2)) 'days)']...
%    , 'VerticalAlignment','Middle','HorizontalAlignment','Center','FontSize',28,'Color','red', 'Interpreter', 'latex')
%    text(18,-7.5,['Prof: ' num2str(round(LZ(nz),2)) ' m']...
%    , 'VerticalAlignment','Middle','HorizontalAlignment','Center','FontSize',28,'Color','red', 'Interpreter', 'latex')
        
    h1=plot([0 60 60 0 0] ,[-20 -20 20 20 -20],'-','color','black', 'LineWidth',1.5)     
%=====================
    hh=colorbar(pl1,'Ticks',[-3:3:3]*10^(-8) ,'TickLabels',...
    {'$-3$' '$0$' '$3$'},'Location','northoutside',...
    'Position',[ 0.07+(0.28+0.02)*(nc-1)+0.15 0.55-(0.4+0.06)*(nf-1)+0.30  0.1  0.01],'TickLabelInterpreter','latex','FontSize',32.0); 
    set(get(hh,'title'),'string',[{'$\mathrm{OW (s^{-2}) \times 10^{-8}}$'}],'FontSize',32.0,'Interpreter','latex','color',rgb('black')); 
    hh.Box= 'on'
    hh.Color= rgb('black')
%==========================FS===============================    
nc=2; nf=1
pos1 = [0.07+(0.28+0.02)*(nc-1) 0.55-(0.4+0.03)*(nf-1) 0.28 0.4];
    pl1=subplot('Position',pos1), hold on    
    pcolor(LXY1/1000, LXY2/1000-100,FF(:,:,nz)), shading interp
    %[hh cc]=contour(LXY1/1000, LXY2/1000,PV(:,:,nz),[-.1,-.1]*10^(-9),'ShowText','off','color',rgb('green'));
    %cc.LineWidth=1.5, cc.LineStyle='-'  ;      
    caxis([-1.0 1.0]*10^(-17));
    colormap(BWR2(360))
    %cmocean('balance','pivot',0)
    %title('Frontogenesis function ($1/s^5$)', 'Interpreter', 'latex','FontSize',bb) %colorbar 
    set(gca,'FontName','Iowan Old Style','FontSize',aa);
    axis equal,xlim([0 60]), ylim([-20 20])
    set(gca, 'box','on','XTickLabel',[ 9 18 27 ],'XTick',[],'YTickLabel',[9 18 27 ]-18,'YTick',[]-18)
h1=plot([0 60 60 0 0] ,[-20 -20 20 20 -20],'-','color','black', 'LineWidth',1.5)        
%=====================
    hh=colorbar(pl1,'Ticks',[-1:1:1]*10^(-17) ,'TickLabels',...
    {'$-1$' '$0$' '$1$'},'Location','northoutside',...
    'Position',[ 0.07+(0.28+0.02)*(nc-1)+0.15 0.55-(0.4+0.06)*(nf-1)+0.30  0.1  0.01],'TickLabelInterpreter','latex','FontSize',32.0); 
    set(get(hh,'title'),'string',[{'$\mathrm{F_s (s^{-5}) \times 10^{-17}}$'}],'FontSize',32.0,'Interpreter','latex','color',rgb('black')); 
    hh.Box= 'on'
    hh.Color= rgb('black')  
                h3=text(1,21,{'$\textbf{b)}$'},'color','black',...
    'FontSize',35.0,'HorizontalAlignment','left','VerticalAlignment','bottom','Interpreter','latex','Rotation',0);
%==== vertical velocity ====================================================    
nc=1; nf=2
pos1 = [0.07+(0.28+0.02)*(nc-1) 0.55-(0.4+0.06)*(nf-1) 0.28 0.4];
    pl1=subplot('Position',pos1), hold on 
    je=pcolor(LXY1/1000, LXY2/1000-100,w(:,:,nz)*60*60*24), shading interp
    je.FaceAlpha=0.7;
    %[hh cc]=contour(LXY1/1000, LXY2/1000,PV(:,:,nz),[-.1,-.1]*10^(-9),'ShowText','off','color',rgb('green'));
    %cc.LineWidth=1.5, cc.LineStyle='-'  
    caxis([-150 150])
    colormap(BWR2(360))
    %cmocean('balance','pivot',0)
    %title('Vertical Velocity ($m/d$)', 'Interpreter','latex','FontSize',bb), %colorbar ;
    set(gca,'FontName','Iowan Old Style','FontSize',aa);
    axis equal,xlim([0 60]), ylim([-20 20])
    set(gca, 'box','on','XTickLabel',[15 30 45 ],'XTick',[15 30 45],'YTickLabel',[-15 0 15 ],'YTick',[-15 0 15]);
    xlabel('X-axis (km)', 'Interpreter', 'latex')%,ylabel('Profundidad (m)', 'Interpreter', 'latex');
    ylabel('Y-axis (km)', 'Interpreter', 'latex')%,ylabel('Profundidad (m)', 'Interpreter', 'latex');
h1=plot([0 60 60 0 0] ,[-20 -20 20 20 -20],'-','color','black', 'LineWidth',1.5)   
            h3=text(1,21,{'$\textbf{d)}$'},'color','black',...
    'FontSize',35.0,'HorizontalAlignment','left','VerticalAlignment','bottom','Interpreter','latex','Rotation',0);
%=====================
    hh=colorbar(pl1,'Ticks',[-150:150:150]*10^(0) ,'TickLabels',...
    {'$-150$' '$0$' '$150$'},'Location','northoutside',...
    'Position',[ 0.07+(0.28+0.02)*(nc-1)+0.15 0.55-(0.4+0.06)*(nf-1)+0.30  0.1  0.01],'TickLabelInterpreter','latex','FontSize',32.0); 
    set(get(hh,'title'),'string',[{'$\mathrm{w (m/d)}$'}],'FontSize',32.0,'Interpreter','latex','color',rgb('black')); 
    hh.Box= 'on'
    hh.Color= rgb('black')     
%=============================================================    
nc=2; nf=2
pos1 = [0.07+(0.28+0.02)*(nc-1) 0.55-(0.4+0.06)*(nf-1) 0.28 0.4];
    pl1=subplot('Position',pos1), hold on 
    pcolor(LXY1/1000, LXY2/1000-100,Tprima(:,:,nz)), shading interp 
    %[hh cc]=contour(LXY1/1000, LXY2/1000,PV(:,:,nz),[-.1,-.1]*10^(-9),'ShowText','off','color',rgb('green'));
    %cc.LineWidth=1.5, cc.LineStyle='-'
    caxis([-0.25 0.25])
    colormap(BWR2(360))
    %cmocean('balance','pivot',0)
%    title('Anomaly Temperature ($^oC$)','Interpreter', 'latex','FontSize',bb);% colorbar 
    set(gca,'FontName','Iowan Old Style','FontSize',aa);
    axis equal,xlim([0 60]), ylim([-20 20])
    set(gca, 'box','on','XTickLabel',[15 30 45],'XTick',[15 30 45],'YTickLabel',[9 18 27 ]-18,'YTick',[]-18)
    xlabel('X-axis (km)', 'Interpreter', 'latex')%,ylabel('Profundidad (m)', 'Interpreter', 'latex');
    h1=plot([0 60 60 0 0] ,[-20 -20 20 20 -20],'-','color','black', 'LineWidth',1.5)  
            h3=text(1,21,{'$\textbf{e)}$'},'color','black',...
    'FontSize',35.0,'HorizontalAlignment','left','VerticalAlignment','bottom','Interpreter','latex','Rotation',0);
%=====================
    hh=colorbar(pl1,'Ticks',[-0.25:0.25:0.25]*10^(0) ,'TickLabels',...
    {'$-0.25$' '$0$' '$0.25$'},'Location','northoutside',...
    'Position',[ 0.07+(0.28+0.02)*(nc-1)+0.15 0.55-(0.4+0.06)*(nf-1)+0.30  0.1  0.01],'TickLabelInterpreter','latex','FontSize',32.0); 
    set(get(hh,'title'),'string',[{'$\mathrm{\theta_a (^0 C)}$'}],'FontSize',32.0,'Interpreter','latex','color',rgb('black')); 
    hh.Box= 'on'
    hh.Color= rgb('white')         
    %hh=quiver(LXY1(1:npx:end,1:npy:end)/1000, LXY2(1:npx:end,1:npy:end)/1000-18,...
    %uu,vv);
    %hh.MarkerSize=5 ;hh.AutoScaleFactor=1.5, hh.MaxHeadSize=1.2;
    %hh.Color='green'%[0.6 0.6 0.6]%[0 50 0]/255%[34 139 34]/255%hh.Color=[47 79 79]/255;
    %hh.LineStyle='-';
%==========================================================================
nc=3; nf=2
pos1 = [0.07+(0.28+0.02)*(nc-1) 0.55-(0.4+0.06)*(nf-1) 0.28 0.4];
    pl1=subplot('Position',pos1), hold on     
    je=pcolor(LXY1/1000, LXY2/1000-100,VHF(:,:,nz)); shading interp
    je.FaceAlpha=0.7;
    caxis([-1000 1000]);
    colormap(BWR2(360));
    %cmocean('balance','pivot',0)
    %title('Vertical Heat Flux ($W/m^2$)', 'Interpreter', 'latex','FontSize',bb), colorbar 
    set(gca,'FontName','Iowan Old Style','FontSize',aa);
    axis equal,xlim([0 60]), ylim([-20 20])
    set(gca, 'box','on','XTickLabel',[15 30 45],'XTick',[15 30 45],'YTickLabel',[9 18 27 ]-18,'YTick',[]-18)
    xlabel('X-axis (km)', 'Interpreter', 'latex')%,ylabel('Profundidad (m)', 'Interpreter', 'latex');
    h1=plot([0 60 60 0 0] ,[-20 -20 20 20 -20],'-','color','black', 'LineWidth',1.5)  
    h3=text(1,21,{'$\textbf{f)}$'},'color','black',...
    'FontSize',35.0,'HorizontalAlignment','left','VerticalAlignment','bottom','Interpreter','latex','Rotation',0);
%=====================
    hh=colorbar(pl1,'Ticks',[-1000:1000:1000]*10^(0) ,'TickLabels',...
    {'$-1000$' '$0$' '$1000$'},'Location','northoutside',...
    'Position',[ 0.07+(0.28+0.02)*(nc-1)+0.15 0.55-(0.4+0.06)*(nf-1)+0.30  0.1  0.01],'TickLabelInterpreter','latex','FontSize',32.0); 
    set(get(hh,'title'),'string',[{'$\mathrm{ VHF (W/m^2)}$'}],'FontSize',32.5,'Interpreter','latex','color',rgb('black')); 
    hh.Box= 'on'
    hh.Color= rgb('black') 
    
%=================================================    
    set(gcf,'Color','white','Renderer','zbuffer')
    set(gca,'Color','white','XColor','black', ...
        'YColor','black')
%-----------------------------------------------------------
end
%==================================================



%=====================================================================================================
%=====================================================================================================
%=====================================================================================================
caso=2; %Forced case with KPP
%=========================
dirdata=['../DATA/EXP';  num2str(caso) '/'];
NT=30:60:(960*60)-30;%NT=30:60:(1440*60)-30;
%=============================================================
close all
f1=figure('units','normalized','Position',[0.3 0.1 0.65 0.65]), hold on
%=====================================================================        
con=0
for nt=NT(700) 
%==============vertical velocity===============
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
%==============anomaly temp===============
[Tprima] = fun_temp_anomaly(dirdata,dirPC,nt,caso);
%=================Okubo-weiss==================
ta=2*(-(ux.*bx + vx.*by)).*bx;
tb=2*(-(uy.*bx + vy.*by)).*by;
FF=ta+tb; % frontogenesis tendency
f0=1e-4;  
%=======PV-ERtel=================================
q1=(wy-vz).*bx;
q2=(uz-wx).*by;
q3=(vx-uy+f0).*bz;
q=q1+q2+q3;
%================================================
Div=ux+vy;  % Divergence
rv=vx-uy;   % Relative Vorticity
ss=vx+uy;   % 
sn=ux-vy;   %
S=sqrt(ss.^2+ sn.^2); % Deformation
OW=sn.^2+ss.^2-rv.^2; %Okubo-weiss
cp=3850;ro=1025;
VHF=cp*ro*w.*Tprima; % vertical heat flux
%================================================
aa=34;bb=26;
%f1=figure('units','normalized','Position',[0.3 0.1 0.65 0.65]), hold on
nz=14; % 50 meters deep.
nc=3; nf=1
pos1 = [0.07+(0.28+0.02)*(nc-1) 0.55-(0.4+0.03)*(nf-1) 0.28 0.4];
    pl1=subplot('Position',pos1)
        hold on
        %============PV========
        pcolor(LXY1/1000, LXY2/1000-100,q(:,:,nz)); shading interp;    
        %[hh cc]=contour(LXY1/1000, LXY2/1000,PV(:,:,nz),[-.1,-.1]*10^(-9),'ShowText','off','color',rgb('green'));
        %cc.LineWidth=1.5, cc.LineStyle='-'
        [hh cc]=contour(LXY1/1000, LXY2/1000-100,q(:,:,nz),[0,0],'ShowText','off','color',[0.6 0.6 0.6]);
        cc.LineWidth=0.8, cc.LineStyle='-' ;  
        caxis([-2 2]*10^(-9));
        colormap(BWR2(360));
        %hh=colorbar
        %cmocean('balance','pivot',0)
    %    title('Ertel-PV ($1/s^3$)', 'Interpreter','latex','FontSize',bb), %colorbar ;
        set(gca,'FontName','Iowan Old Style','FontSize',aa);
        axis equal,xlim([0 60]), ylim([-20 20])
        set(gca, 'box','on','XTickLabel',[ 9 18 27 ],'XTick',[],'YTickLabel',[9 18 27 ]-18,'YTick',[]-18)
        h3=text(1,21,{'$\textbf{c)}$'},'color','black',...
    'FontSize',35.0,'HorizontalAlignment','left','VerticalAlignment','bottom','Interpreter','latex','Rotation',0);
%=====================
    h1=plot([0 60 60 0 0] ,[-20 -20 20 20 -20],'-','color','black', 'LineWidth',1.5)       
%=====================
    hh=colorbar(pl1,'Ticks',[-2:2:2]*10^(-9) ,'TickLabels',...
    {'$-2$' '$0$' '$2$'},'Location','northoutside',...
    'Position',[ 0.07+(0.28+0.02)*(nc-1)+0.15 0.55-(0.4+0.06)*(nf-1)+0.30  0.1  0.01],'TickLabelInterpreter','latex','FontSize',32.0); 
    set(get(hh,'title'),'string',[{'$\mathrm{q (s^{-3}) \times 10^{-9}}$'}],'FontSize',32.0,'Interpreter','latex','color',rgb('black')); 
    hh.Box= 'on'
    hh.Color= rgb('black')  
    colormap(BWR2(360));
%==========================================================================    
nc=1; nf=1
pos1 = [0.07+(0.28+0.02)*(nc-1) 0.55-(0.4+0.03)*(nf-1) 0.28 0.4];
    pl1=subplot('Position',pos1), hold on
    pcolor(LXY1/1000, LXY2/1000-100,OW(:,:,nz)), shading interp
    %[hh cc]=contour(LXY1/1000, LXY2/1000,PV(:,:,nz),[-.1,-.1]*10^(-9),'ShowText','off','color',rgb('green'));
    %cc.LineWidth=1.5, cc.LineStyle='-'
    caxis([-3 3]*10^(-8))
    colormap(BWR2(360))
    %cmocean('balance','pivot',0)
    %title('Okubo-Weiss ($1/s^2$)', 'Interpreter', 'latex','FontSize',bb)% colorbar 
    set(gca,'FontName','Iowan Old Style','FontSize',aa);
        axis equal,xlim([0 60]), ylim([-20 20])
    set(gca, 'box','on','XTickLabel',[ 9 18 27 ],'XTick',[],'YTickLabel',[-15 0 15],'YTick',[-15 0 15])
    ylabel('Y-axis (km)', 'Interpreter', 'latex')%,ylabel('Profundidad (m)', 'Interpreter', 'latex');
    dt=30
            h3=text(1,21,{'$\textbf{a)}$'},'color','black',...
    'FontSize',35.0,'HorizontalAlignment','left','VerticalAlignment','bottom','Interpreter','latex','Rotation',0);
%    text(18,-5,['Time: ' num2str(round(nt*dt/(60*60),2)) ' hours (' num2str(round(nt*dt/(60*60*24),2)) 'days)']...
%    , 'VerticalAlignment','Middle','HorizontalAlignment','Center','FontSize',28,'Color','red', 'Interpreter', 'latex')
%    text(18,-7.5,['Prof: ' num2str(round(LZ(nz),2)) ' m']...
%    , 'VerticalAlignment','Middle','HorizontalAlignment','Center','FontSize',28,'Color','red', 'Interpreter', 'latex')
        
    h1=plot([0 60 60 0 0] ,[-20 -20 20 20 -20],'-','color','black', 'LineWidth',1.5)     
%=====================
    hh=colorbar(pl1,'Ticks',[-3:3:3]*10^(-8) ,'TickLabels',...
    {'$-3$' '$0$' '$3$'},'Location','northoutside',...
    'Position',[ 0.07+(0.28+0.02)*(nc-1)+0.15 0.55-(0.4+0.06)*(nf-1)+0.30  0.1  0.01],'TickLabelInterpreter','latex','FontSize',32.0); 
    set(get(hh,'title'),'string',[{'$\mathrm{OW (s^{-2}) \times 10^{-8}}$'}],'FontSize',32.0,'Interpreter','latex','color',rgb('black')); 
    hh.Box= 'on'
    hh.Color= rgb('black')
%==========================FS===============================    
nc=2; nf=1
pos1 = [0.07+(0.28+0.02)*(nc-1) 0.55-(0.4+0.03)*(nf-1) 0.28 0.4];
    pl1=subplot('Position',pos1), hold on    
    pcolor(LXY1/1000, LXY2/1000-100,FF(:,:,nz)), shading interp
    %[hh cc]=contour(LXY1/1000, LXY2/1000,PV(:,:,nz),[-.1,-.1]*10^(-9),'ShowText','off','color',rgb('green'));
    %cc.LineWidth=1.5, cc.LineStyle='-'  ;      
    caxis([-1.0 1.0]*10^(-17));
    colormap(BWR2(360))
    %cmocean('balance','pivot',0)
    %title('Frontogenesis function ($1/s^5$)', 'Interpreter', 'latex','FontSize',bb) %colorbar 
    set(gca,'FontName','Iowan Old Style','FontSize',aa);
    axis equal,xlim([0 60]), ylim([-20 20])
    set(gca, 'box','on','XTickLabel',[ 9 18 27 ],'XTick',[],'YTickLabel',[9 18 27 ]-18,'YTick',[]-18)
h1=plot([0 60 60 0 0] ,[-20 -20 20 20 -20],'-','color','black', 'LineWidth',1.5)        
%=====================
    hh=colorbar(pl1,'Ticks',[-1:1:1]*10^(-17) ,'TickLabels',...
    {'$-1$' '$0$' '$1$'},'Location','northoutside',...
    'Position',[ 0.07+(0.28+0.02)*(nc-1)+0.15 0.55-(0.4+0.06)*(nf-1)+0.30  0.1  0.01],'TickLabelInterpreter','latex','FontSize',32.0); 
    set(get(hh,'title'),'string',[{'$\mathrm{F_s (s^{-5}) \times 10^{-17}}$'}],'FontSize',32.0,'Interpreter','latex','color',rgb('black')); 
    hh.Box= 'on'
    hh.Color= rgb('black')  
                h3=text(1,21,{'$\textbf{b)}$'},'color','black',...
    'FontSize',35.0,'HorizontalAlignment','left','VerticalAlignment','bottom','Interpreter','latex','Rotation',0);
%==== vertical velocity ====================================================    
nc=1; nf=2
pos1 = [0.07+(0.28+0.02)*(nc-1) 0.55-(0.4+0.06)*(nf-1) 0.28 0.4];
    pl1=subplot('Position',pos1), hold on 
    je=pcolor(LXY1/1000, LXY2/1000-100,w(:,:,nz)*60*60*24), shading interp
    je.FaceAlpha=0.7;
    %[hh cc]=contour(LXY1/1000, LXY2/1000,PV(:,:,nz),[-.1,-.1]*10^(-9),'ShowText','off','color',rgb('green'));
    %cc.LineWidth=1.5, cc.LineStyle='-'  
    caxis([-150 150])
    colormap(BWR2(360))
    %cmocean('balance','pivot',0)
    %title('Vertical Velocity ($m/d$)', 'Interpreter','latex','FontSize',bb), %colorbar ;
    set(gca,'FontName','Iowan Old Style','FontSize',aa);
    axis equal,xlim([0 60]), ylim([-20 20])
    set(gca, 'box','on','XTickLabel',[15 30 45 ],'XTick',[15 30 45],'YTickLabel',[-15 0 15 ],'YTick',[-15 0 15]);
    xlabel('X-axis (km)', 'Interpreter', 'latex')%,ylabel('Profundidad (m)', 'Interpreter', 'latex');
    ylabel('Y-axis (km)', 'Interpreter', 'latex')%,ylabel('Profundidad (m)', 'Interpreter', 'latex');
h1=plot([0 60 60 0 0] ,[-20 -20 20 20 -20],'-','color','black', 'LineWidth',1.5)   
            h3=text(1,21,{'$\textbf{d)}$'},'color','black',...
    'FontSize',35.0,'HorizontalAlignment','left','VerticalAlignment','bottom','Interpreter','latex','Rotation',0);
%=====================
    hh=colorbar(pl1,'Ticks',[-150:150:150]*10^(0) ,'TickLabels',...
    {'$-150$' '$0$' '$150$'},'Location','northoutside',...
    'Position',[ 0.07+(0.28+0.02)*(nc-1)+0.15 0.55-(0.4+0.06)*(nf-1)+0.30  0.1  0.01],'TickLabelInterpreter','latex','FontSize',32.0); 
    set(get(hh,'title'),'string',[{'$\mathrm{w (m/d)}$'}],'FontSize',32.0,'Interpreter','latex','color',rgb('black')); 
    hh.Box= 'on'
    hh.Color= rgb('black')     
%=============================================================    
nc=2; nf=2
pos1 = [0.07+(0.28+0.02)*(nc-1) 0.55-(0.4+0.06)*(nf-1) 0.28 0.4];
    pl1=subplot('Position',pos1), hold on 
    pcolor(LXY1/1000, LXY2/1000-100,Tprima(:,:,nz)), shading interp 
    %[hh cc]=contour(LXY1/1000, LXY2/1000,PV(:,:,nz),[-.1,-.1]*10^(-9),'ShowText','off','color',rgb('green'));
    %cc.LineWidth=1.5, cc.LineStyle='-'
    caxis([-0.25 0.25])
    colormap(BWR2(360))
    %cmocean('balance','pivot',0)
%    title('Anomaly Temperature ($^oC$)','Interpreter', 'latex','FontSize',bb);% colorbar 
    set(gca,'FontName','Iowan Old Style','FontSize',aa);
    axis equal,xlim([0 60]), ylim([-20 20])
    set(gca, 'box','on','XTickLabel',[15 30 45],'XTick',[15 30 45],'YTickLabel',[9 18 27 ]-18,'YTick',[]-18)
    xlabel('X-axis (km)', 'Interpreter', 'latex')%,ylabel('Profundidad (m)', 'Interpreter', 'latex');
    h1=plot([0 60 60 0 0] ,[-20 -20 20 20 -20],'-','color','black', 'LineWidth',1.5)  
            h3=text(1,21,{'$\textbf{e)}$'},'color','black',...
    'FontSize',35.0,'HorizontalAlignment','left','VerticalAlignment','bottom','Interpreter','latex','Rotation',0);
%=====================
    hh=colorbar(pl1,'Ticks',[-0.25:0.25:0.25]*10^(0) ,'TickLabels',...
    {'$-0.25$' '$0$' '$0.25$'},'Location','northoutside',...
    'Position',[ 0.07+(0.28+0.02)*(nc-1)+0.15 0.55-(0.4+0.06)*(nf-1)+0.30  0.1  0.01],'TickLabelInterpreter','latex','FontSize',32.0); 
    set(get(hh,'title'),'string',[{'$\mathrm{\theta_a (^0 C)}$'}],'FontSize',32.0,'Interpreter','latex','color',rgb('black')); 
    hh.Box= 'on'
    hh.Color= rgb('white')         
    %hh=quiver(LXY1(1:npx:end,1:npy:end)/1000, LXY2(1:npx:end,1:npy:end)/1000-18,...
    %uu,vv);
    %hh.MarkerSize=5 ;hh.AutoScaleFactor=1.5, hh.MaxHeadSize=1.2;
    %hh.Color='green'%[0.6 0.6 0.6]%[0 50 0]/255%[34 139 34]/255%hh.Color=[47 79 79]/255;
    %hh.LineStyle='-';
%==========================================================================
nc=3; nf=2
pos1 = [0.07+(0.28+0.02)*(nc-1) 0.55-(0.4+0.06)*(nf-1) 0.28 0.4];
    pl1=subplot('Position',pos1), hold on     
    je=pcolor(LXY1/1000, LXY2/1000-100,VHF(:,:,nz)); shading interp
    je.FaceAlpha=0.7;
    caxis([-1000 1000]);
    colormap(BWR2(360));
    %cmocean('balance','pivot',0)
    %title('Vertical Heat Flux ($W/m^2$)', 'Interpreter', 'latex','FontSize',bb), colorbar 
    set(gca,'FontName','Iowan Old Style','FontSize',aa);
    axis equal,xlim([0 60]), ylim([-20 20])
    set(gca, 'box','on','XTickLabel',[15 30 45],'XTick',[15 30 45],'YTickLabel',[9 18 27 ]-18,'YTick',[]-18)
    xlabel('X-axis (km)', 'Interpreter', 'latex')%,ylabel('Profundidad (m)', 'Interpreter', 'latex');
    h1=plot([0 60 60 0 0] ,[-20 -20 20 20 -20],'-','color','black', 'LineWidth',1.5)  
    h3=text(1,21,{'$\textbf{f)}$'},'color','black',...
    'FontSize',35.0,'HorizontalAlignment','left','VerticalAlignment','bottom','Interpreter','latex','Rotation',0);
%=====================
    hh=colorbar(pl1,'Ticks',[-1000:1000:1000]*10^(0) ,'TickLabels',...
    {'$-1000$' '$0$' '$1000$'},'Location','northoutside',...
    'Position',[ 0.07+(0.28+0.02)*(nc-1)+0.15 0.55-(0.4+0.06)*(nf-1)+0.30  0.1  0.01],'TickLabelInterpreter','latex','FontSize',32.0); 
    set(get(hh,'title'),'string',[{'$\mathrm{ VHF (W/m^2)}$'}],'FontSize',32.5,'Interpreter','latex','color',rgb('black')); 
    hh.Box= 'on'
    hh.Color= rgb('black') 
    
%=================================================    
    set(gcf,'Color','white','Renderer','zbuffer')
    set(gca,'Color','white','XColor','black', ...
        'YColor','black')
%-----------------------------------------------------------
end
%==================================================




