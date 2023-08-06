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
#=====extra analysis====================================
    mpl.rcParams['axes.linewidth'] = 2.5
    mpl.rc('xtick',labelsize=20)
    mpl.rc('ytick',labelsize=20)
    mpl.rc('text',usetex=False)
    mpl.rcParams['xtick.direction'] = 'out'
    mpl.rcParams['ytick.direction'] = 'out'
    mpl.rcParams['font.family']='Times New Roman'
#==================
# For Unforced case ==================================
    test=1
#====================    
    f  = 1./17.4
    if test==1:name='coVHFtest1r.npz'
    if test==4:name='coVHFtest4r.npz'
    datos=np.load(dirdata+name)
    Eiso=datos['Eiso'];Eiso=Eiso[:,:-1];
    om=datos['om'];om=om[:-1];kiso=datos['kiso']
    EmVHF=Eiso.T*kiso[None,...]*om[...,None];
    kVHFkpp=(EmVHF).sum(axis=0) 
    fVHFkpp=(EmVHF).sum(axis=1)     
    #========================== 
    if test==1: name='coVHFtest5r.npz'
    if test==4: name='coVHFtest6r.npz'
    datos=np.load(dirdata+name)
    Eiso=datos['Eiso'];Eiso=Eiso[:,:-1];
    om=datos['om'];om=om[:-1];kiso=datos['kiso']
    EmVHF=Eiso.T*kiso[None,...]*om[...,None];
    kVHF=(EmVHF).sum(axis=0) 
    fVHF=(EmVHF).sum(axis=1)  
#===========================
    plt.close('all')
    fig = plt.figure(figsize=(13,6.5))
    ax1 = plt.subplot2grid((10,15),(1,1),rowspan=8, colspan=6);       
    ax1.plot(kiso,cp * ro *kVHFkpp,'-',linewidth=2.5,color='black',label='KPP active')
    ax1.plot(kiso,cp * ro *kVHF,'--',linewidth=2.5,color='black',label='KPP inactive')
    ax1.set_xscale('log')
    ax1.set_title('(a)',horizontalalignment='right',size=24)  
    ax1.set_ylabel('$\sum_{\omega} \kappa \omega VHF(\kappa,\omega) [W/m^2]$',size=20)
    ax1.set_xlabel('$\kappa$: Horizontal Wavenumber [cpkm]',size=22)
    ax1.plot([LS,LS],[-1000,2000],'--',linewidth=1.4,color='gray')
    if test==1: ax1.plot([1/3,1/3],[-200,200],'-',linewidth=1.6,color='green')
    #if test==4: ax1.plot([1/5,1/5],[-200,200],'-',linewidth=1.6,color='green')
    if test==4:ax1.plot([1/2,1/2],[-200,200],'-',linewidth=1.6,color='green')
    ax1.plot([1/35,2.5],[0,0],'--',linewidth=1.4,color='gray')    
    ax1.set_ylim(-600,1400)
    ax1.set_xlim(1/35.,2.5)
    ax1a=ax1.twiny()
    ax1a.set_xscale('log')
    ax1a.set_xticks([1./30,1./10,LS,1])
    ax1a.set_xticklabels(['30','10','$\mathrm{L_s}$','1'])
    ax1a.set_xlim(1/35.,2.5)
    ax1a.set_xlabel('Wavelength [km]',size=22)
    #===========================
    ax1 = plt.subplot2grid((10,15),(1,10),rowspan=8, colspan=4);       
    ax1.plot(cp * ro *fVHFkpp,om,'-',linewidth=2.5,color='black',label='KPP active')
    ax1.plot(cp * ro *fVHF,om,'--',linewidth=2.5,color='black',label='KPP inactive')
    ax1.set_yscale('log')
    ax1.set_title('(b)',horizontalalignment='right',size=24)  
    ax1.set_xlabel('$\sum_{\kappa} \kappa \omega VHF(\kappa,\omega) [W/m^2]$',size=20)
    ax1.set_ylabel('$\omega$: Frequency [cph]',size=22)
    ax1.legend(fontsize=15)
    ax1.plot([0,0],[1/(10*24),1/1],'--',linewidth=1.4,color='gray')
    ax1.set_xlim(-200,700)
    if test==1: ax1.plot([-100,100],[1/10,1/10],'-',linewidth=1.6,color='green')
    if test==4: ax1.plot([-100,100],[1/6,1/6],'-',linewidth=1.6,color='green')
    ax1a=ax1.twiny()
    ax1a.set_xticks([0])
    ax1a.set_xticklabels(['1'],color='white')
    ax1a.set_xlabel('--',size=16,color='white')
    #ax1.legend(fontsize=14)
    ax1.set_ylim(1/(10*24),1/1)
    ax1a=ax1.twinx()
    ax1a.set_yscale('log')
    ax1a.set_yticks([1/(7*24),1/24,f,1/6])
    ax1a.set_yticklabels(['7','1','$\mathrm{f}$','6h'])
    ax1a.set_ylim(1/(10*24),1/1)
    ax1a.set_ylabel('Period [days]',size=22,rotation=-90,
        verticalalignment='bottom')
    #============================================  



#==================
# To Forced case ==================================
    test=4
#====================  
    f  = 1./17.4
    if test==1:name='coVHFtest1r.npz'
    if test==4:name='coVHFtest4r.npz'
    datos=np.load(dirdata+name)
    Eiso=datos['Eiso'];Eiso=Eiso[:,:-1];
    om=datos['om'];om=om[:-1];kiso=datos['kiso']
    EmVHF=Eiso.T*kiso[None,...]*om[...,None];
    kVHFkpp=(EmVHF).sum(axis=0) 
    fVHFkpp=(EmVHF).sum(axis=1)     
    #========================== 
    if test==1: name='coVHFtest5r.npz'
    if test==4: name='coVHFtest6r.npz'
    datos=np.load(dirdata+name)
    Eiso=datos['Eiso'];Eiso=Eiso[:,:-1];
    om=datos['om'];om=om[:-1];kiso=datos['kiso']
    EmVHF=Eiso.T*kiso[None,...]*om[...,None];
    kVHF=(EmVHF).sum(axis=0) 
    fVHF=(EmVHF).sum(axis=1)  
#===========================
    plt.close('all')
    fig = plt.figure(figsize=(13,6.5))
    ax1 = plt.subplot2grid((10,15),(1,1),rowspan=8, colspan=6);       
    ax1.plot(kiso,cp * ro *kVHFkpp,'-',linewidth=2.5,color='black',label='KPP active')
    ax1.plot(kiso,cp * ro *kVHF,'--',linewidth=2.5,color='black',label='KPP inactive')
    ax1.set_xscale('log')
    ax1.set_title('(a)',horizontalalignment='right',size=24)  
    ax1.set_ylabel('$\sum_{\omega} \kappa \omega VHF(\kappa,\omega) [W/m^2]$',size=20)
    ax1.set_xlabel('$\kappa$: Horizontal Wavenumber [cpkm]',size=22)
    ax1.plot([LS,LS],[-1000,2000],'--',linewidth=1.4,color='gray')
    if test==1: ax1.plot([1/3,1/3],[-200,200],'-',linewidth=1.6,color='green')
    #if test==4: ax1.plot([1/5,1/5],[-200,200],'-',linewidth=1.6,color='green')
    if test==4:ax1.plot([1/2,1/2],[-200,200],'-',linewidth=1.6,color='green')
    ax1.plot([1/35,2.5],[0,0],'--',linewidth=1.4,color='gray')    
    ax1.set_ylim(-600,1400)
    ax1.set_xlim(1/35.,2.5)
    ax1a=ax1.twiny()
    ax1a.set_xscale('log')
    ax1a.set_xticks([1./30,1./10,LS,1])
    ax1a.set_xticklabels(['30','10','$\mathrm{L_s}$','1'])
    ax1a.set_xlim(1/35.,2.5)
    ax1a.set_xlabel('Wavelength [km]',size=22)
    #===========================
    ax1 = plt.subplot2grid((10,15),(1,10),rowspan=8, colspan=4);       
    ax1.plot(cp * ro *fVHFkpp,om,'-',linewidth=2.5,color='black',label='KPP active')
    ax1.plot(cp * ro *fVHF,om,'--',linewidth=2.5,color='black',label='KPP inactive')
    ax1.set_yscale('log')
    ax1.set_title('(b)',horizontalalignment='right',size=24)  
    ax1.set_xlabel('$\sum_{\kappa} \kappa \omega VHF(\kappa,\omega) [W/m^2]$',size=20)
    ax1.set_ylabel('$\omega$: Frequency [cph]',size=22)
    ax1.legend(fontsize=15)
    ax1.plot([0,0],[1/(10*24),1/1],'--',linewidth=1.4,color='gray')
    ax1.set_xlim(-200,700)
    if test==1: ax1.plot([-100,100],[1/10,1/10],'-',linewidth=1.6,color='green')
    if test==4: ax1.plot([-100,100],[1/6,1/6],'-',linewidth=1.6,color='green')
    ax1a=ax1.twiny()
    ax1a.set_xticks([0])
    ax1a.set_xticklabels(['1'],color='white')
    ax1a.set_xlabel('--',size=16,color='white')
    #ax1.legend(fontsize=14)
    ax1.set_ylim(1/(10*24),1/1)
    ax1a=ax1.twinx()
    ax1a.set_yscale('log')
    ax1a.set_yticks([1/(7*24),1/24,f,1/6])
    ax1a.set_yticklabels(['7','1','$\mathrm{f}$','6h'])
    ax1a.set_ylim(1/(10*24),1/1)
    ax1a.set_ylabel('Period [days]',size=22,rotation=-90,
        verticalalignment='bottom')
    #============================================   