%%==============Fig2================================
clear all, close all
load ('../DATA/ParamSeries.mat')
f1=figure('units','normalized','Position',[0.3 0.1 0.6 0.5]), hold on
%===================================================
%===================================================
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
%===================================================
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

%===================================================
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
%===================================================
%===========
%===================================================
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
%===================================================
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
%===================================================
%===================================================
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
%===================================================
set(gcf,'Color','white','Renderer','zbuffer');
set(gca,'Color','white','XColor','black', ...
    'YColor','black');
%===================================================
%name=['Fig2'];
%set(f1,'Units','Inches');
%pos = get(f1,'Position');
%set(f1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
%print(f1,name,'-dpdf','-r350')
%print(f1,name,'-dpng','-r350')
%% ================================end ==========================