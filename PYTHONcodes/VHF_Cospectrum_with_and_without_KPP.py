    import numpy as np
    from numpy import  pi
    import os
    import warnings
    import glob
    import scipy.io as sci
    import matplotlib as mpl
    #mpl.use('Agg')
    from matplotlib.patches import Patch
    import matplotlib.pyplot as plt
    import matplotlib.mlab as mlab
    import matplotlib.cm as cm
    from matplotlib import colors, ticker, cm
    from matplotlib.colors import LinearSegmentedColormap
    import wf_spectrum
    #import spectral_kinematics_flux_heat
    import matplotlib.mathtext as mathtext
    from matplotlib.colors import LogNorm
    from scipy import signal
    import h5py
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    from math import log
    #===============================================
    import h5py
    import mat73
    import SPECTRAL as SP    
    from scipy.io import loadmat
    #===============================================

    dirdata='../DATA/'

    #=================VHF=====================================
    # DATA explain
    # coVHFtest1r.npz   Unforced VHF Co-spectrum with KPP  
    # coVHFtest4r.npz   Forced VHF Co-spectrum with KPP
    # coVHFtest5r.npz   Unforced VHF Co-spectrum without KPP
    # coVHFtest6r.npz   Forced VHF Co-spectrum without KPP
    #================
    # VWtest1r.npz      Unforced vertical velocity spectrum with KPP          
    # VWtest4r.npz      Forced vertical velocity spectrum with KPP
    # VWtest5r.npz      Unforced vertical velocity spectrum without KPP
    # VWtest6r.npz      Forced vertical velocity spectrum without KPP
    #=======================================================
    normwt = mpl.colors.BoundaryNorm(bounds, cmap.N)
    cp=3850
    ro=1025
    LS=1/2.3
    mpl.rcParams['axes.linewidth'] = 2.5
    mpl.rc('xtick',labelsize=18)
    mpl.rc('ytick',labelsize=18)
    mpl.rc('text',usetex=False)
    mpl.rcParams['xtick.direction'] = 'out'
    mpl.rcParams['ytick.direction'] = 'out'
    mpl.rcParams['font.family']='Times New Roman'
    
    plt.close('all')
    
    #fig = plt.figure(figsize=(9,6.0))
    # esce=1 to plot with KPP
    # esce=2 to plot without KPP
    for esce in range(1,3):
        fig = plt.figure(figsize=(9,6.0))
        for nf in range(1,3):
        print (nf)
        cmap=plt.cm.get_cmap('seismic')
    
        if nf==1:
            ax1 = ax1 = plt.subplot2grid((12,17),(1,1),rowspan=8, colspan=7);       
            if esce==1:
                name='coVHFtest1r.npz';name2='VWtest1r.npz';cmap=plt.cm.get_cmap('seismic')  ; fun='VHF'  
            if esce==2:
               name='coVHFtest5r.npz';name2='VWtest5r.npz';cmap=plt.cm.get_cmap('seismic')  ; fun='VHF'  
        if nf==2:
            ax1  = plt.subplot2grid((12,17),(1,9),rowspan=8, colspan=7)
            if esce==1:
                name='coVHFtest4r.npz';name2='VWtest4r.npz';cmap=plt.cm.get_cmap('seismic'); fun='VHF'
            if esce==2:
               name='coVHFtest6r.npz';name2='VWtest6r.npz';cmap=plt.cm.get_cmap('seismic')  ; fun='VHF'  
        #==================Co-VHF===================
        datos=np.load(dirdata+name)
        Eiso=datos['Eiso'];Eiso=Eiso[:,:-1];
        om=datos['om'];om=om[:-1];kiso=datos['kiso']
        #EmVHF=Eiso.T; 
        EmVHF=Eiso.T*kiso[None,...]*om[...,None];
        cmap=plt.cm.get_cmap('seismic')
        cs = ax1.pcolormesh(kiso, om, cp * ro * EmVHF, cmap=cmap, shading='flat',vmin=-50, vmax=50)
        #===================Vel. Vertical==========
        datos=np.load(dirdata+name2)
        Eiso=datos['Eiso'];Eiso=Eiso[:,:-1];
        om=datos['om'];om=om[:-1];
        #Em=Eiso.T;     
        Em=Eiso.T*kiso[None,...]*om[...,None];
        valmax=Em.max()
        # ====  Figure =====
        #cs = ax1.pcolormesh(kiso, om, cp * ro * EmVHF, cmap=cmap, shading='flat')#,vmin=-100, vmax=100)
        levels = [0.01,0.05, 0.1 ,0.2,0.3, 0.4, 0.6, 0.8, 1]
        cs2 = ax1.contour(kiso, om, Em/valmax,levels=levels,cmap='YlGn',vmin=-0.3, vmax=0.75,linewidths=0.7)
             
        
        ax1.set_yscale('log')
        ax1.set_xscale('log')
        ax1.set_ylim(1/(24*10.),1.)
        ax1.set_xlim(1/35.,2.5)
        ax1a=ax1.twiny()
        ax1a.set_yscale('log')
        ax1a.set_xscale('log')
        ax1a.set_xticks([1./30,1./10,1/5.,LS,1])
        ax1a.set_xticklabels(['30','10','','$\mathrm{L_s}$','1'],size=18)  
        ax1a.set_xlim(1/35.,2.5)
        ax1a=ax1.twinx()
        ax1a.set_yscale('log')
        ax1a.set_yticks([1./(24*7),1./48,1/24.,1/12])
        ax1a.set_yticklabels(['','','',''])
        ax1a.set_ylim(1/(24*10.),1.)
        if nf==2:
            ax1a.set_yticks([1./(24*7),1./48,1/24.])
            ax1a.set_yticklabels(['7','2','1'],size=18)
            ax1.set_yticks([]) ; ax1.set_yticklabels([])
        ks = np.array([1.e-3,2.5])
        ks = np.array([1.e-3,2.5])
        f  = 1./17.4
        D1= 1/24.0
        ax1.plot(ks,[f,f],'k--',linewidth=1.2,color='gray')
        #ax1.text(1.5,f+0.0025,r'$f$',color='green',size=20)
        ax1.plot(ks,[D1,D1],'k--',linewidth=1.2,color='black')
        ax1.plot([LS,LS],[1/(24*10),1],'k--',linewidth=1.2,color='black') 
        #ax1.text(1/15,1/(2.2),(letra[nf-1]+')'),color='black',size=14)
        if nf==1: ax1.text(1/30,1/(1.83),'a)',color='black' ,size=26)
        if nf==2: ax1.text(1/30,1/(1.83),'b)',color='black' ,size=26)
        
        if nf==1: ax1.text(1/150,1/(18),'$\mathrm{\omega}$: Frequency [cph]',color='black',size=22,rotation=90,verticalalignment='center')
        if nf==2: ax1.text(1/0.25,1/(18),"Period [days]",color='black',size=22,rotation=-90,verticalalignment='center')   
        if nf==2: ax1.text(1/900,1/(34*24),'$\mathrm{\kappa}$: Horizontal wavenumber [cpkm]',color='black',size=22)    
        if nf==2: ax1.text(1/300,1/(0.015*24),"Wavelength [km]",color='black',size=22)
            #cbar = fig.colorbar(cs, ax=ax1, orientation='horizontal', pad=0.2)
            #cbar.set_label('cp * ro * EmVHF')
        if nf==2:                                  
            cbar_ax = fig.add_axes([0.645,0.07,0.3,0.03])
            aaa=fig.colorbar(cs, cax=cbar_ax,
                         norm=mpl.colors.Normalize(),
                         ticks=[-80,-40,0,40,80],orientation = 'horizontal')
            aaa.set_ticklabels(['','-40','0','40',''])
            ax1.text(1/0.85,1/(24*68),r'$\mathrm{\kappa \omega }$ Re $[\mathrm{\rho C_p}$ $\mathrm{\widehat{w}\widehat{\theta_a}]}}$ $~\mathrm{(W/m^2)}$',color='black',size=16,horizontalalignment='center')
    