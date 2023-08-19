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
#================================
#DIVE1r.npz   Unforced divergence
#VORT1r.npz   Unforced vorticity
#DIVE4r.npz   Forced divergence
#VORT4r.npz   Forced vorticity
#==================
# dimension for the graphics #===========================
nnx=54; nte=1/5; nte2=1/40
ccx=0.7 ; ccy=0.1 ; lltx=1/7 ; llty=1/(10*90); sizele=18;
#=======================================================
LS=1/2.3
letra='abcdefghijklmn'
#============================
mpl.rcParams['axes.linewidth'] = 2
mpl.rc('xtick',labelsize=18)
mpl.rc('ytick',labelsize=18)
mpl.rc('text',usetex=False)
mpl.rcParams['xtick.direction'] = 'out'
mpl.rcParams['ytick.direction'] = 'out'
mpl.rcParams['font.family']='Times New Roman'
#============================================================


plt.close('all')
fig = plt.figure(figsize=(11,6.0))
#============================================================#============================================================
nny=9 ;nnx=46
for nf in [1,2,3,4]:#range(1,15):
    print (nf)
    cmap=plt.cm.get_cmap('ocean_r')
    #cmap=plt.cm.get_cmap('gist_ncar_r')
    if nf==1:
        ax1 = plt.subplot2grid((nny,nnx),(1,0),rowspan=6, colspan=10);
        name='DIVE1r.npz'  ; fun='$\mathrm{\hat{\delta}}$' ; let='a)'        
    if nf==2:
        ax1 = plt.subplot2grid((nny,nnx),(1,11),rowspan=6, colspan=10);        
        name='VORT1r.npz' ; fun='$\mathrm{\hat{\zeta}}$'  ; let='b)'       
    if nf==3:
        ax1 = plt.subplot2grid((nny,nnx),(1,25),rowspan=6, colspan=10);        
        name='DIVE4r.npz'; fun='$\mathrm{\hat{\delta}}$' ; let='c)'   
    if nf==4:
        ax1 = plt.subplot2grid((nny,nnx),(1,36),rowspan=6, colspan=10);        
        name='VORT4r.npz'; fun='$\mathrm{\hat{\zeta}}$'  ; let='d)' 

    
    #======================================================
    # to get the maximun value (Forced case) to normalize
    datos=np.load(dirdata+'VORT4r.npz')
    Eiso=datos['Eiso'];
    Eiso=Eiso[:,:-1];
    om=datos['om'];om=om[:-1];
    kiso=datos['kiso'];
    Em=Eiso.T*om[...,None];
    valmax=Em.max()       
    #=========================================    
    datos=np.load(dirdata+name)
    Eiso=datos['Eiso'];Eiso=Eiso[:,:-1];
    Em=Eiso.T*om[...,None]; 
    if nf==13 or nf==14: 
        print('::::::::: Hello :::::::::')
    else:Em=Em/valmax
    #=========================================
    cs=plt.pcolormesh(kiso,om,np.log(Em), cmap=cmap,shading='flat')
    ax1.set_yscale('log')
    ax1.set_xscale('log')
    ax1.set_ylim(1/(24*10.),1.)
    ax1.set_xlim(1/30.,2.5)
    ax1a=ax1.twiny()
    ax1a.set_yscale('log')
    ax1a.set_xscale('log')
    if nf==1 or nf==2 or nf==3 or nf==4 or nf==9 or nf==11 or nf==13:
        ax1a.set_xticks([1/30,1./18,1./10,1/5.,LS,1])
        ax1a.set_xticklabels(['30','','10','','$\mathrm{L_s}$','1'])    
    ax1a.set_xlim(1/30.,2.5)
    ax1a=ax1.twinx()
    ax1a.set_yscale('log')
    ax1a.set_yticks([1./(24*5),1./48,1/24.,1/12])
    ax1a.set_yticklabels(['','','',''])
    ax1a.set_ylim(1/(24*10.),1.)
    if nf==2 or nf==4:
        ax1a.set_yticks([1./(24*7),1./48,1/24.])
        ax1a.set_yticklabels(['7','2','1'],size=18)
    ks = np.array([1.e-3,2.5])
    ks = np.array([1.e-3,2.5])
    f  = 1./17.4
    D1= 1/24.0
    if nf==13 or nf ==14:
        corl='black'; corl2='black'
    else:
        corl='red'; corl2='magenta'
    ax1.text(1/20,1/(1.9),let,color='black',size=22)    
    ax1.plot(ks,[f,f],'k--',linewidth=0.7,color='white')
    ax1.plot(ks,[1/(24),1/(24)],'k--',linewidth=1.2,color='red')
    #ax1.text(1.5,f+0.0025,r'$f$',color=corl,size=sizele)
        #ax1.plot(ks,[D1,D1],'k--',linewidth=1.2,color='forestgreen')
    ax1.plot([LS,LS],[1/(24*10),1],'k--',linewidth=1.2,color='red')    
    ax1.text(1/0.9,1/(6.5*24),fun,color='black',size=25)
   
    #ax1.text(1/15,1/(2.2),(letra[nf-1]+')'),color='black',size=14)
    
    if nf==5 or nf==13 or nf==5 or nf==7 or nf==9 or nf==11 or nf==13: 
        ax1.set_xticks([]) ; ax1.set_xticklabels([])
    if nf==3 or nf==4 or nf==2 or nf==6 or nf==7 or nf==8 or nf==9 or nf==10:
        ax1.set_yticks([]) ; ax1.set_yticklabels([])
    if nf==11 or nf==12 or nf==13 or nf==14: ax1.set_yticks([]) ; ax1.set_yticklabels([])    
#============================================================#============================================================
    if nf==1:ax1.set_ylabel("$\mathrm{\omega}$: Frequency [cph]",color='black',size=19,rotation=90,verticalalignment='center')
    if nf==3: 
        ax1.text(1/60,1/(0.31),"Wavelength [km]",color='black',size=19,horizontalalignment='center')
    if nf==3: ax1.text(nte2,1/(37*24),'$\mathrm{\kappa}$: Horizontal wavenumber [cpkm]',color='black',size=19,horizontalalignment='center')
    if nf==4: ax1.text(1/0.20,1/(17),"Period [days]",color='black',size=18,rotation=-90,verticalalignment='center')
#============================================================#============================================================
    if nf<13: clim = [-8,0]; plt.clim(clim)     
    if nf==13 or nf==14: clim = [-60,60]; plt.clim(clim)
    if nf==4: 
        scale=10**(0)   
        cbar_ax = fig.add_axes([0.735,ccy,0.2,0.02])
        aaa=fig.colorbar(cs, cax=cbar_ax,
                     norm=mpl.colors.Normalize(),
                     ticks=[-6,-3,0],orientation = 'horizontal')
        aaa.set_ticklabels(['$\mathrm{10^{-6}}$','$\mathrm{10^{-3}}$','$\mathrm{1}$'])
        ax1.text(1/2,1/(40*25),r'$\mathrm{\kappa \omega|\widehat{C}|^2}$/max($\mathrm{\kappa \omega |\widehat{C}|^2}$)',color='black',size=15,horizontalalignment='center')
#===============================Save grafico !!!!!!!!!!!!=======================================