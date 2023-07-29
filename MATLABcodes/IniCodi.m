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
nt=30;
dirdata=['../DATA/']
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
%%================================================
% Done Figure
%%================================================
con2=1
close all
f1=figure('units','normalized','Position',[0.1 0.1 0.5 0.7]); % ===========
pos1 = [0.7 0.75 0.25 0.15];
    pl1=subplot('Position',pos1)
    Usup(:)=U(1,:,1);
    plot(LY/1000-100,Usup,'-','color',rgb('Darkgray'),'linewidth',3.5)
    xlim([-30 30]); ylim([-0.01 0.25]);     set(gca,'FontName','Iowan Old Style','FontSize',ziseleter);
    xlabel({'Y-axis (km)'},'FontSize',ziseleter,'Interpreter', 'latex','Rotation',0);
    ylabel({'$u$ ($m/s)$'},'FontSize',ziseleter,'Interpreter', 'latex');
    set(gca,'XTick',[-15 0 15],'XTickLabel',{'-15' '0' '15'})
    set(gca,'YTick',[0.0 0.1 0.15 0.2],'YTickLabel',{'0' '0.1' '' '0.2'}),
    h3=text(22.5,0.174,{'$\textbf{b)}$'},'color','black',...
    'FontSize',35.0,'HorizontalAlignment','left','VerticalAlignment','bottom','Interpreter','latex','Rotation',0);
%=====================================================================
pos1 = [0.642 0.15 0.16 0.45];
    pl1=subplot('Position',pos1)
    roo(:,1,1)=nanmean(nanmean(RHO,1),2);
    plot(roo-1000,LZ,'-','color',rgb('Darkgray'),'linewidth',3.5)
    xlim([25.45 26.5]);
    ylim([-290 -1.9]);     set(gca,'FontName','Iowan Old Style','FontSize',ziseleter);
    ylabel('Depth(m)','FontSize',26, 'Interpreter', 'latex','Rotation',90);
    xlabel({'$\sigma_{\theta} -10^3$ ($Kg/m^{3}$)'},'FontSize',28,'Interpreter', 'latex');
    set(gca,'YTick',[-250 -200 -150 -100 -50 ],'ZTickLabel',{'-250' '-200' '-150' '-100' '-50'});
    h3=text(26.3,-32,{'$\textbf{c)}$'},'color','black',...
    'FontSize',35.0,'HorizontalAlignment','left','VerticalAlignment','bottom','Interpreter','latex','Rotation',0);
%=====================================================================
pos1 = [0.83 0.15 0.16 0.45];
    pl1=subplot('Position',pos1)
    N2(:,1,1)=nanmean(nanmean(bz,1),2);
    plot(N2*10^5,LZ,'-','color',rgb('Darkgray'),'linewidth',3.5)
    xlim([-0.5 5.5]); ylim([-290 -1.9]);     set(gca,'FontName','Iowan Old Style','FontSize',ziseleter);
    %ylabel('Depth(m)','FontSize',26, 'Interpreter', 'latex','Rotation',90);
    xlabel({'$N^2$ ($s^{-2})$ $\times 10^{-5}$'},'FontSize',26,'Interpreter', 'latex');
    set(gca,'YTick',[ ],'YTickLabel',[ ])
    set(gca,'YTick',[-250 -200 -150 -100 -50 -2],'ZTickLabel',{'' ' ' ' ' ' ' ' ' ' '});
    h3=text(4.3,-32,{'$\textbf{d)}$'},'color','black',...
    'FontSize',35.0,'HorizontalAlignment','left','VerticalAlignment','bottom','Interpreter','latex','Rotation',0);
%=====================================================================
set(gcf,'Color','white','Renderer','zbuffer');
set(gca,'Color','white','XColor','black', ...
    'YColor','black');
%===================================================================== 
pos1 = [0.1+(con2-1)*0.40 0.09 0.4 0.85];
    pl1=subplot('Position',pos1)
        hold on
        h=slice(Y3-100,X3,Z3,THETA,[-30 30],[0 59.4],[-1.9]); shading interp % aqui    
        %h(3).FaceAlpha=0.6;
        %h2(nl).FaceAlpha=1;
        %h2(nl).FaceAlpha=1;
        set(gca,'Ydir','reverse')%,'YTick',[9 18 27],'YTickLabel',[-9 0 9]);
        set(gca,'DataAspectRatio',[1 1 2.5])
        view([-29 21]);
        colormap(jet(300))
        caxis([14.5 15.25])
        zlim([-250 -1.9]),xlim([-30 30]),ylim( [0 60])
        grid off
        box_obj = gca;
        set(box_obj, 'Box', 'on', 'BoxStyle', 'full', 'LineWidth', 1.1);
        
        
     L=0.45
    hh=colorbar(pl1,'Ticks',[14.50 14.75 15.00 15.25] ,'TickLabels',...
    {'$14.50$' '$14.75$' '$15.00$' '$15.25$'},'Location','eastoutside',...
    'Position',[ 0.510 (1-L)/2+0.19  0.035/3.5 L],'TickLabelInterpreter','latex','FontSize',32.0); 
    set(get(hh,'title'),'string',[{'$^0 C$'}],'FontSize',32.0,'Interpreter','latex'); 
    cbarrow down       
%================== Temperature contour
        lvls=[0.04:0.06:0.25]
        h2=contourslice(Y3-100,X3,Z3,U,[-30 30],[0 59.4],[-1.9],lvls);
        for nl=1:length(h2)
            h2(nl).EdgeColor=rgb('white');
            h2(nl).FaceAlpha=1;
            h2(nl).LineWidth=3.0;
        end        
%================== boyance contour
        lvls=[0.0211:0.001:0.02905]
        h2=contourslice(Y3-100,X3,Z3,b,[-30 30],[0 59.4],[-1.9],lvls);
        for nl=1:length(h2)
            h2(nl).EdgeColor=rgb('Darkgray');
            h2(nl).FaceAlpha=1;
            h2(nl).LineWidth=1.5;
            h2(nl).LineStyle='--'
        end       
%================== boyance contour
        lvls=[0.02911 0.02944 0.02977]
        h2=contourslice(Y3-100,X3,Z3,b,[-30 30],[0 59.4],[-1.9],lvls);
        for nl=1:length(h2)
            h2(nl).EdgeColor=rgb('black');
            h2(nl).FaceAlpha=1;
            h2(nl).LineWidth=1.5;
            h2(nl).LineStyle='--'
        end           
%==================General aspects       
        set(gca,'FontName','Iowan Old Style','FontSize',ziseleter);
        q=ylabel({'X-axis','(km)'},'FontSize',ziseleter,'Interpreter', 'latex','Rotation',-40);
        q.HorizontalAlignment='right'
        q=xlabel({'Y-axis','(km)'},'FontSize',ziseleter,'Interpreter', 'latex','Rotation',14);
        q.HorizontalAlignment='left'
        zlabel('Depth(m)','FontSize',ziseleter, 'Interpreter', 'latex','Rotation',90);
        set(gca,'YTick',[15 30 45],'YTickLabel',{'15' '30' '45'})
        set(gca,'XTick',[-15 0 15],'XTickLabel',{'-15' '0' '15'})
        set(gca,'ZTick',[-200 -150 -100 -50 -2],'ZTickLabel',{'-200' '-150' '-100' '-50' ''});
        h3=text(-18,0,20,{'$\textbf{a)}$'},'color','black',...
    'FontSize',35.0,'HorizontalAlignment','left','VerticalAlignment','bottom','Interpreter','latex','Rotation',0);
%=====================================================================