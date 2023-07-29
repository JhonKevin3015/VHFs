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
%% ====VHF in 3D ====================
%================================================
close all
f1=figure('units','normalized','Position',[0.1 0.1 0.62 0.74]); 
%===================================

con2=0
for caso=[3 4]; 
dirdata=['../DATA/EXP';  num2str(caso) '/'];
nt=NT(700);
dias=nt*dt/(60*60*24);
CAMPO='W';   
dat=rdmds([dirdata CAMPO],nt);
dat(:,1:2,:)=NaN;
dat(:,end-1:end,:)=NaN;
w=dat;
%===================================
u=rdmds([dirdata 'U'],nt);
u(:,1:2,:)=NaN;
u(:,end-1:end,:)=NaN;
[uy ux uz]=gradient(u,LY,LX,LZ);
%===================================
v=rdmds([dirdata 'V'],nt);
v(:,1:2,:)=NaN;
v(:,end-1:end,:)=NaN;
[vy vx vz]=gradient(v,LY,LX,LZ);
%========================================= 
%==============anomaly temp===============
[Tprima] = fun_temp_anomaly(dirdata,dirPC,nt,caso);
cp=3850;ro=1025;
VHF=cp*ro*w.*Tprima;
%=========================================
[by,bx,bz,b] = fun_bonyance_and_gradient(dirdata,nt,LX,LY,LZ);
RI=bz./(uz.*uz+vz.*vz);
g=9.81;
rho0=1028.60;
RHO=-1028.60*b/g+rho0;
%=====================================================================
con2=con2+1;
%===================================================================== 
%close all
%f1=figure('units','normalized','Position',[0.1 0.1 0.6 0.74]); % ==========

   pos1 = [0.07+(con2-1)*0.38 0.11 0.38 0.78];
    pl1=subplot('Position',pos1)
        hold on
    %h=slice(Y3-18,X3,Z3,THETA,[-17.6 17],[0 35.4],[-1.9]); shading interp % aqui   
    h=slice(Y3-100,X3,Z3,VHF,[],[],[-2.0 -50 -100 -150 -200 -250]); shading interp % aqui    
    h2=slice(Y3-100,X3,Z3,VHF,[9] ,[0],[]); shading interp % aqui 
        %h(3).FaceAlpha=0.6;
        %h2(nl).FaceAlpha=1;
        %h2(nl).FaceAlpha=1;
        set(gca,'Ydir','reverse')%,'YTick',[9 18 27],'YTickLabel',[-9 0 9]);
        set(gca,'DataAspectRatio',[1 1 2.5])
        view([-31 18]);
        caxis([-1000 1000])
        colormap(BWR2(360))
        zlim([-200 -2.0]),xlim([-20 20]),ylim( [0 60])
        box on, grid off
        ax = gca;
        ax.BoxStyle = 'full'; 
        ax.LineWidth= 1.7
h(1).FaceAlpha=0.6;
h(2).FaceAlpha=0.8;
h(3).FaceAlpha=0.8;
h(4).FaceAlpha=0.8;
h(5).FaceAlpha=0.8;
h(6).FaceAlpha=0.8;
h2(1).FaceAlpha=0.4;
h2(2).FaceAlpha=0.4;
            
%================== boyance contour
        lvls=[0.0211:0.0005:0.028]
        h2=contourslice(Y3-100,X3,Z3,b,[] ,[0],[],lvls);
        for nl=1:length(h2)
            h2(nl).EdgeColor=rgb('green');
            h2(nl).FaceAlpha=1;
            h2(nl).LineWidth=1.0;
            h2(nl).LineStyle='--'
        end       
% %================== boyance contour
        lvls=[0.028 :0.0001: 0.03]
        h2=contourslice(Y3-100,X3,Z3,b,[] ,[0],[-2 -50 -100 -150 -200],lvls);
        for nl=1:length(h2)
            h2(nl).EdgeColor=rgb('green');
            h2(nl).FaceAlpha=0.6;
            h2(nl).LineWidth=1;
            h2(nl).LineStyle='-'
        end           
%==================General aspects       
        ziseleter=30;
        set(gca,'FontName','Iowan Old Style','FontSize',ziseleter);
        q=ylabel({'X-axis','(km)'},'FontSize',ziseleter,'Interpreter', 'latex','Rotation',-40);
        q.HorizontalAlignment='right'
        q=xlabel({'Y-axis','(km)'},'FontSize',ziseleter,'Interpreter', 'latex','Rotation',14);
        q.HorizontalAlignment='left'
        zlabel('Depth(m)','FontSize',ziseleter, 'Interpreter', 'latex','Rotation',90);
        set(gca,'YTick',[15 30 45 ],'YTickLabel',{'15' '30' '45'})
        set(gca,'XTick',[-15 0 15],'XTickLabel',{'-15' '0' '15'})
if con2==1;        h3=text(-30,5,5,{'$\textbf{a)}$'},'color','black',...
    'FontSize',40.0,'HorizontalAlignment','left','VerticalAlignment','bottom','Interpreter','latex','Rotation',0);
    set(gca,'ZTick',[-200 -150 -100 -50 -2],'ZTickLabel',{'-200' '-150' '-100' '-50' ''});
    zlabel('Depth(m)','FontSize',ziseleter, 'Interpreter', 'latex','Rotation',90);
    VHFper1(:)=nanmean(reshape(VHF(:,401:600,:),[300*200 60]),1);    
%h3=text(-18,40,-225,{'$\textbf{Unforced}$'},'color','black',...
    %'FontSize',40.0,'HorizontalAlignment','center','VerticalAlignment','bottom','Interpreter','latex','Rotation',0);
    BZper1(:)=nanmean(reshape(bz(:,401:600,:),[300*200 60]),1);
    RIper1(:)=nanmean(reshape(RI(:,401:600,:),[300*200 60]),1);
end
if con2==2;         h3=text(-30,5,5,{'$\textbf{b)}$'},'color','black',...
    'FontSize',40.0,'HorizontalAlignment','left','VerticalAlignment','bottom','Interpreter','latex','Rotation',0);
    set(gca,'ZTick',[-200 -150 -100 -50 -2],'ZTickLabel',[]);
    zlabel([],'FontSize',ziseleter, 'Interpreter', 'latex','Rotation',90);
    VHFper2(:)=nanmean(reshape(VHF(:,401:600,:),[300*200 60]),1);
%h3=text(-18,40,-225,{'$\textbf{Forced}$'},'color','black',...
    %'FontSize',40.0,'HorizontalAlignment','center','VerticalAlignment','bottom','Interpreter','latex','Rotation',0);       
    BZper2(:)=nanmean(reshape(bz(:,401:600,:),[300*200 60]),1);
    RIper2(:)=nanmean(reshape(RI(:,401:600,:),[300*200 60]),1);
end

set(gcf,'Color','white','Renderer','zbuffer');
set(gca,'Color','white','XColor','black', ...
    'YColor','black');
end

    hh=colorbar(pl1,'Ticks',[-1000:500:1000] ,'TickLabels',...
    {'$-1000$' '$-500$' '$0$' '$500$' '$1000$'},'Location','northoutside',...
    'Position',[0.33  0.89  0.25  0.01],'TickLabelInterpreter','latex','FontSize',32.0); 
    set(get(hh,'title'),'string',[{'$\textbf{VHF} \mathrm{(W/m^2)}$'}],'FontSize',32.0,'Interpreter','latex'); 
%=====================================================================

   %figure () 
   corr=rgb('green');%[0.4 0.4 0.4]
   pos1 = [0.86 0.16 0.13 0.605];
   pl1=subplot('Position',pos1), hold on
   p1=plot(VHFper1,LZ,'--','color',rgb('black'),'LineWidth',3.5)
   p2=plot(VHFper2,LZ,'-','color','black','LineWidth',3.5)
   p3=plot(BZper1*10^6, LZ,'--','color',corr,'LineWidth',2.5)
   p4=plot(BZper2*10^6, LZ,'-','color',corr,'LineWidth',2.5)
   plot([-15 50 50 -15 -15],[-200 -200 -2 -2 -200],'-','color','black','LineWidth',2.1)
   plot(LZ*0,LZ,'-','color','black','LineWidth',1)
   set(gca,'ZTick',[-200 -150 -100 -50 -2],'ZTickLabel',{'-200' '-150' '-100' '-50' ''});
   set(gca,'FontName','Iowan Old Style','FontSize',ziseleter);
   ylim([-200 -2]), xlim([-15 50])
   box on, grid off
   xlabel([{'$\overline{\textbf{VHF}}$' ,'$\mathrm{(W/m^2)}$'}],'FontSize',ziseleter, 'Interpreter', 'latex','Rotation',0);
   h3=text(-25,0,{'$\textbf{c)}$'},'color','black',...
    'FontSize',40.0,'HorizontalAlignment','left','VerticalAlignment','bottom','Interpreter','latex','Rotation',0);
   text(0,0,{'${0}$'},'color',corr,...
    'FontSize',30.0,'HorizontalAlignment','center','VerticalAlignment','bottom','Interpreter','latex','Rotation',0);
   text(20,0,{'${2}$'},'color',corr,...
    'FontSize',30.0,'HorizontalAlignment','center','VerticalAlignment','bottom','Interpreter','latex','Rotation',0);
   text(40,0,{'${4}$'},'color',corr,...
    'FontSize',30.0,'HorizontalAlignment','center','VerticalAlignment','bottom','Interpreter','latex','Rotation',0);
   text(15,11,{'$\overline{b_z}$ (10$^{-5}$/s$^2$)'},'color',corr,...
    'FontSize',30.0,'HorizontalAlignment','center','VerticalAlignment','bottom','Interpreter','latex','Rotation',0);
%         [hLg, icons]=legend([p1 p2],'$\textbf{Forced}$','$\textbf{Unforced}$','Interpreter', 'latex'...
%         ,'Orientation','vertical','FontSize',30);
%         icons = findobj(icons,'Type','patch');
%         icons = findobj(icons,'Marker','none','-xor');
%         hLg.Position=[0.89 0.81 0.08 0.06],hLg.Box='off'
%         hLg.FontSize=23; 
%===================================       
        hh=legend('Unforced','Forced','Unforced','Forced');
        hh.FontSize=22
        hh.Position=[0.86 0.9 0.08 0.06]; 
        hh.Box='off'
%=====================================================