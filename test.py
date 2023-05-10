#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 12 16:26:57 2021


@author: jhon
"""
import numpy as np
import matplotlib.pyplot as plt


dx=0.2
dy=0.2
xs=np.arange(-15,15,dx) # km
ys=np.arange(-15,15,dy)
xm,ym=np.meshgrid(xs,ys) 
k=2*pi/5
l=2*pi/5


t=np.arange(0,20,1/24) # 20*24 pasos de timpo% 20 dias 
T1=1
w1=2*pi/T1

u = np.zeros((len(x1),len(x1),len(t)))
v = np.zeros((len(x1),len(x1),len(t)))

for nt in np.arange(0,len(t),1):
    print(nt)
    Z1 = np.sin(k*xm-w1*t[nt])
    Z2 = np.sin(l*ym-w1*t[nt])
    u[:,:,nt]=Z1[:]    
    v[:,:,nt]=Z2[:] 
#    print(np)


aa=u[0,0,:]
bb=v[0,0,:]
   
%matplotlib qt
plt.figure()
plt.plot(t,aa)
plt.plot(t,bb)
plt.show()
    

%matplotlib qt
fig = plt.figure(figsize= (9,5))
ax = fig.gca(projection='3d')
surface=ax.plot_surface(xm,ym,Z1,cmap="hot")
fig.colorbar(surface)
plt.show()


fig1=plt.figure (1)
ax1=plt.subplot(2,1,1)
p1=ax1.pcolor(xm,ym,Z1,cmap="hot")
ax2=plt.subplot(2,1,2)
p2=ax2.pcolor(xm,ym,Z2,cmap="hot")
fig1.colorbar(p2)
plt.show()

def update_plot(frame_number, ZT3d, plot):
    plot[0].remove()
    plot[0] = ax.plot_surface(xm, ym, ZT3d[:,:,frame_number], cmap="magma")

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

fps = 10 # frame per sec
frn = 480 # frame number of the animation
plot = [ax.plot_surface(xm, ym, ZT3d[:,:,0], color='0.75', rstride=1, cstride=1)]
#ax.set_zlim(0,1.1)
ani = animation.FuncAnimation(fig, update_plot, frn, fargs=(ZT3d, plot), interval=10/fps)

fn = 'plot_surface_animation_funcanimation'
ani.save(fn+'.mp4',writer='ffmpeg',fps=fps)
ani.save(fn+'.gif',writer='imagemagick',fps=fps)
#==========================END of animations====================
#=======================TO co-stectrum  ========================
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
import matplotlib.mathtext as mathtext
from matplotlib.colors import LogNorm
from scipy import signal
#========================================
import h5py
import mat73
import SPECTRAL as SP
#=======================================================
# cambio entre ODL-vertion and 73-vertion 
from scipy.io import loadmat
dirdata='/Users/oceano/DATA-to-pythonOLD/'
#dirdata='/Users/oceano/DATA-to-python/'
#name1='tp-exp1-2.mat'; name2='wp-exp1-2.mat'
name1='tp-exp1.mat'; name2='wp-exp1.mat'
#===================
#del tp; del wt; del WT
anot = loadmat(dirdata+name1)
#anot = mat73.loadmat(dirdata+name1)
tp=anot['Tp'][:,]
anot = loadmat(dirdata+name2) 
#anot = mat73.loadmat(dirdata+name2)
wp=anot['Wp'][:,] 
#    WT=wp*tp   
#===================    
#tp = signal.detrend(tp,axis=0,type='linear')
#tp = signal.detrend(tp,axis=1,type='linear')
#tp = signal.detrend(tp,axis=2,type='linear')
#wp = signal.detrend(wp,axis=0,type='linear')
#wp = signal.detrend(wp,axis=1,type='linear')
#wp = signal.detrend(wp,axis=2,type='linear')
print('------- End detrending --------')
#=============================
cp=3850
ro=1025
[a,b,c]=wp.shape
x1=np.int64(a/2)
x2=np.int64(b/2)
plt.figure(), plt.plot(tp[x1,x2,:])
plt.figure(), plt.plot(wp[x1,x2,:]*24*60*60)
plt.figure(), plt.plot(cp*ro*wp[x1,x2,:]*tp[x1,x2,:])

LX = np.linspace(0,a*0.2, a)
LY = np.linspace(0,b*0.2, b)
LT= np.linspace(0,c/48, c)
nti=913-1
wp2d=wp[:,:,nti]
tp2d=tp[:,:,nti]


plt.close('all')
fig = plt.figure(13,figsize=(12,4))
ax3=plt.subplot(133)
cs=plt.pcolor(LX,LY,ro*cp*wp2d.T*tp2d.T,cmap=plt.cm.get_cmap('seismic'),shading='flat')
clim = [-1000,1000];plt.clim(clim);plt.colorbar(cs,label=r'wp*tp*Cp*ro(W/m2)')
ax3.set_aspect('equal', 'box')

ax1=plt.subplot(131)
cs=plt.pcolor(LX,LY,wp2d.T*24*60*60,cmap=plt.cm.get_cmap('seismic'),shading='flat')
clim = [-100,100];plt.clim(clim);plt.colorbar(cs,label=r'');ax1.set_title ('V. velocity`')
ax1.set_aspect('equal', 'box')

ax2=plt.subplot(132)
cs=plt.pcolor(LX,LY,tp2d.T,cmap=plt.cm.get_cmap('seismic'),shading='flat')
plt.colorbar(cs,label=r''), ax2.set_title ('A. Temp')
ax2.set_aspect('equal', 'box');clim = [-0.3,0.3];plt.clim(clim)
#========
[a,b,c]=wp.shape
#===================================================
#del coE
d1=0.2;d2=0.2;d3=0.5
con=0
wind='hann'#'kais'#'kais'#'flat'#'hann'# 'flat'
for ncor in range(1,11):
    con=con+1
    print (ncor)
    pdt=48
    ventana=10 # dias    
    aa=np.arange(0+pdt*(ncor-1),ventana*24*2+pdt*(ncor-1))
    wpc=wp[:,:,aa] 
    tpc=tp[:,:,aa]
#    wpc = signal.detrend(wpc,axis=2,type='linear')
#    tpc = signal.detrend(tpc,axis=2,type='linear')
    if ncor == 1:
        coE,f1,f2,f3,df1,df2,df3=SP.cospec_ab_funtion(wpc,tpc,d1,d2,d3,wind)
    else:        
        coEa,f1,f2,f3,df1,df2,df3=SP.cospec_ab_funtion(wpc,tpc,d1,d2,d3,wind)
        coE=coE+coEa
coE=coE/con
#===================================================
f3.shape
f2.shape
#===================================================
k=f1;l=f2;om=f3
#Eiso = np.empty((212,om.size)) # el 53 te lo sugiere % 212
Eiso = np.empty((64,om.size)) # el 53 te lo sugiere % 212
# save the 2D spectrum
#np.savez(fname,k=k,l=l,Eu=Eu)
for i in range(om.size-1):
    kiso,Eiso[:,i] = SP.calc_ispec(k,l,coE[:,:,i])
#===================================================
np.savez(dirdata+'wt'+name1[6:7]+'r.npz',Eiso=Eiso,kiso=kiso,om=om)
#np.savez(dirdata+'wt'+name1[7:12]+'r.npz',Eiso=Eiso,kiso=kiso,om=om)
print('------- End save --------')
print('------- Finish --------')
#=========================================
omega=om[:-1]
kiso=kiso
E1=Eiso[:,:-1]
cp=3850
ro=1025
LS=1/2.3
#=========================================
import cmocean
cmap = cmocean.cm.balance
#=========================================
bounds1 = np.linspace(-2, 0, 128)
bounds2 = np.linspace(0, 80,128)
bounds=np.concatenate((bounds1, bounds2), axis=None)
#=========================================
norm1 = mpl.colors.BoundaryNorm(bounds, cmap.N)
#=========================================
bounds1 = np.linspace(-80, 0, 128)
bounds2 = np.linspace(0, 80,128)
bounds=np.concatenate((bounds1, bounds2), axis=None)
#=========================================
norm2 = mpl.colors.BoundaryNorm(bounds, cmap.N)
#=========================================

ax2 = plt.subplot2grid((4,13),(0,6),rowspan=2, colspan=6)
ax3 = plt.subplot2grid((4,13),(2,0),rowspan=2, colspan=6)
ax4 = plt.subplot2grid((4,13),(2,6),rowspan=2, colspan=6)

cmap=plt.cm.get_cmap('seismic')

import matplotlib
matplotlib.rcParams['axes.linewidth'] = 2
matplotlib.rc('xtick',labelsize=14)
matplotlib.rc('ytick',labelsize=14)
matplotlib.rc('text',usetex=False)
matplotlib.rcParams['xtick.direction'] = 'out'
matplotlib.rcParams['ytick.direction'] = 'out'
matplotlib.rcParams['font.family']='Times New Roman'

#=================================================
#plt.close('all')
fig = plt.figure(figsize=(6,6))
ax1 = plt.subplot2grid((4,13),(0,0),rowspan=2, colspan=6)
cs=plt.pcolormesh(kiso,om[:-1],cp*ro*E1.T, cmap=cmap,shading='flat')  
#cs=plt.pcolor(kiso,omega,
#              (cp*ro*E1.T),
#              cmap=cmap,shading='flat',norm=norm2)
#              cmap=cmap,shading='flat')  # 'gouraud'

ax1.set_yscale('log')
ax1.set_xscale('log')
ax1.set_ylim(1/(24*5.),1.)
ax1.set_xlim(1/60.,2.5)
ax1.set_yticks([])
ax1.set_yticklabels([])
ax1a=ax1.twiny()
ax1a.set_yscale('log')
ax1a.set_xscale('log')
ax1a.set_xticks([1./50,1./10,1/5,LS,1,2])
ax1a.set_xticklabels(['50','10','5','Ls','1','0.5'])
ax1a.set_yticklabels([])
ax1a.set_xlim(1/60.,2.5)
ax1a=ax1.twinx()
ax1a.set_yscale('log')
ax1a.set_yticks([1./(24*5),1./48,1/24.,1/12])
ax1a.set_yticklabels(['5','2','1','.5'])
#ax1a.set_ylabel('Period [days]',size=18)
ax1a.set_ylim(1/(24*5.),1.)
plt.show()
#ax1a.set_xticklabels([])
#### Some relevant frequencies
ks = np.array([1.e-3,2.5])
#f = coriolis(20.411) ## < ========== Latitude from the name
f  = 1./17.4
K1= 1/24.0
#ax1.plot(ks,[f,f],'k--',linewidth=1.2,color='green')
#ax1.text(1.5,f+0.0025,r'$f$',color='green',size=23)
ax1.plot(ks,[K1,K1],'k--',linewidth=1.2,color='black')
ax1.plot([LS,LS],[1/(24*10),1],'k--',linewidth=1.2,color='w')

ax1.text(1.,1/(24*7.5),r'$Base-expt$',color='black',size=15)


cbar_ax = fig.add_axes([0.91,.11,0.01, 0.8])
fig.colorbar(cs, cax=cbar_ax,
           label=r'coespectro-wp-tp')

clim = [-120,120]
plt.clim(clim)

#===============================       
dirsave='/Users/oceano/MIT-JOB-TESIS/FIGURES/'
plt.savefig(dirsave+'testpy.png', dpi=500, facecolor='w', edgecolor='w',
        orientation='portrait', papertype=None, format='png',
        transparent=False, bbox_inches=None, pad_inches=0.1,
        frameon=None, metadata=None)
#============================================================== 
dirdata='/Users/oceano/DATA-to-pythonOLD/'
name='st-exp1.npz'
name=datos=np.load(dirdata+name)
Eiso=datos['Eiso'];Eiso=Eiso[:,:-1];
kiso=datos['kiso'];
om=datos['om'];om=om[:-1];
#==============================================================
#=================VHF=====================================
cmap = cmocean.cm.balance
bounds1 = np.linspace(-1.0*100, 0, 128)
bounds2 = np.linspace(0, 1.0*100,128)
bounds=np.concatenate((bounds1, bounds2), axis=None)
#=======================================================
normwt = mpl.colors.BoundaryNorm(bounds, cmap.N)
#=======================================================
mpl.rcParams['axes.linewidth'] = 2.5
mpl.rc('xtick',labelsize=18)
mpl.rc('ytick',labelsize=18)
mpl.rc('text',usetex=False)
mpl.rcParams['xtick.direction'] = 'out'
mpl.rcParams['ytick.direction'] = 'out'
mpl.rcParams['font.family']='Times New Roman'



plt.close('all')
fig = plt.figure(figsize=(9,5.4))

for nf in range(1,3):
    print (nf)
    cmap=plt.cm.get_cmap('seismic')

    if nf==1:
        ax1 = ax1 = plt.subplot2grid((12,17),(1,1),rowspan=8, colspan=7);       
        name='wt1r.npz';norm=normwt;cmap=plt.cm.get_cmap('seismic')  ; fun='VHF'  
    if nf==2:
        ax1  = plt.subplot2grid((12,17),(1,9),rowspan=8, colspan=7)
        name='wt3r.npz';norm=normwt;cmap=plt.cm.get_cmap('seismic'); fun='VHF'

    datos=np.load(dirdata+name)
    Eiso=datos['Eiso'];Eiso=Eiso[:,:-1];
    EmVHF=Eiso.T; 
    #cs=plt.pcolormesh(kiso,om,cp*ro*Em, cmap=cmap,shading='gouraud',norm=norm) 
    cs=plt.pcolormesh(kiso,om,cp*ro*EmVHF, cmap=cmap,shading='flat',norm=norm) 
    cs2 = ax.contour(kiso, om, Em/valmax, cmap='viridis',vmin=0, vmax=1)         
    ax1.set_yscale('log')
    ax1.set_xscale('log')
    ax1.set_ylim(1/(24*5.),1.)
    ax1.set_xlim(1/18.,2.5)
    ax1a=ax1.twiny()
    ax1a.set_yscale('log')
    ax1a.set_xscale('log')
    ax1a.set_xticks([1./18,1./10,1/5.,LS,1])
    ax1a.set_xticklabels(['','10','','$L_s$','1'],size=18)  
    ax1a.set_xlim(1/18.,2.5)
    ax1a=ax1.twinx()
    ax1a.set_yscale('log')
    ax1a.set_yticks([1./(24*5),1./48,1/24.,1/12])
    ax1a.set_yticklabels(['','','',''])
    ax1a.set_ylim(1/(24*5.),1.)
    if nf==2:
        ax1a.set_yticks([1./(24*5),1./48,1/24.])
        ax1a.set_yticklabels(['5','2','1'],size=18)
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
    if nf==1: ax1.text(1/16,1/(1.7),'a)',color='black' ,size=24)
    if nf==2: ax1.text(1/16,1/(1.7),'b)',color='black' ,size=24)
    
    if nf==1: ax1.text(1/70,1/(2.5*24),'$\omega$: Frequency [cph]',color='black',size=22,rotation=90)
    if nf==2: ax1.text(1/0.25,1/(2*24),"Period [days]",color='black',size=22,rotation=-90)    
    if nf==2: ax1.text(1/370,1/(15*24),'$\kappa$: Horizontal wavenumber [cpkm]',color='black',size=22)    
    if nf==2: 
        ax1.text(1/100,1/(0.018*24),"Wavelength [km]",color='black',size=22)
        cbar_ax = fig.add_axes([0.645,0.07,0.3,0.03])
        aaa=fig.colorbar(cs, cax=cbar_ax,
                     norm=mpl.colors.Normalize(),
                     ticks=[-100,-50,0,50,100],orientation = 'horizontal')
        aaa.set_ticklabels(['-100','-50','0','50','100'])
        ax1.text(1/0.85,1/(24*30),r'$\kappa$ x $\omega$ x $Re$[$\rho$ $C_p$ $\widehat{w}\widehat{\theta_a}]$ $~(W/m^2$)',color='black',size=16,horizontalalignment='center')

    
dirsave='/Users/oceano/MIT-JOB-TESIS/FIGURES/'
plt.savefig(dirsave+'Fig7paper1CoSpect.png', dpi=800, facecolor='w', edgecolor='w',
        orientation='portrait', papertype=None, format='png',
        transparent=False, bbox_inches=None, pad_inches=0.1,
        frameon=None, metadata=None)
#========================================================================

plt.close('all')
fig = plt.figure(figsize=(9,6.4))

for nf in range(1,3):
    print (nf)
    cmap=plt.cm.get_cmap('seismic')

    if nf==1:
        ax1 = plt.subplot2grid((12,17),(1,1),rowspan=8, colspan=7);       
        name='wt1r.npz';norm=normwt;cmap=plt.cm.get_cmap('seismic')  ; fun='VHF'  
    if nf==2:
        ax1  = plt.subplot2grid((12,17),(1,9),rowspan=8, colspan=7)
        name='wt3r.npz';norm=normwt;cmap=plt.cm.get_cmap('seismic'); fun='VHF'

    datos=np.load(dirdata+name)
    Eiso=datos['Eiso'];Eiso=Eiso[:,:-1];
    EmVHF=Eiso.T;

    if nf==1:        
        name='Vw-exp1.npz';norm=normN  ; fun='w' 
    if nf==2:       
        name='Vw-exp3.npz';norm=normN; fun='w'

    datos=np.load(dirdata+name)
    #datos=np.load(dirdata+name[0:6]+'3'+name[7:])
    Eiso=datos['Eiso'];Eiso=Eiso[:,:-1];
    Em=Eiso.T*kiso[None,...]*om[...,None];
    valmax=Em.max()
    #valmax=5.487623590973079e-05;


    cmap=plt.cm.get_cmap('seismic')

    #fig, ax = plt.subplots()
    # Create the contour plot with all levels
    if nf==1: 
        levels = [0.01, 0.1 ,0.2,0.3, 0.4, 0.6, 0.8, 1]
        cs = ax1.contour(kiso, om, Em/valmax,levels=levels,cmap='YlGn',vmin=-0.3, vmax=1.8,linewidths=0.7)   
    if nf==2: 
        levels = [0.01,0.05, 0.1 ,0.2,0.3, 0.4, 0.6, 0.8, 1]
        cs = ax1.contour(kiso, om, Em/valmax,levels=levels,cmap='YlGn',vmin=-0.3, vmax=0.75,linewidths=0.7)
    
    # Create the pseudocolor plot
    cs2 = ax1.pcolormesh(kiso, om, cp * ro * EmVHF, cmap=cmap, shading='flat',vmin=-100, vmax=100)
    
    # Set the y-axis and x-axis scales to logarithmic scales
    ax1.set_yscale('log')
    ax1.set_xscale('log')
    
    ax1.set_ylim(1/(24*5.),1.)
    ax1.set_xlim(1/18.,2.5)
    
    # Add a colorbar
    #fig.colorbar(cs2)
    
    ax1a=ax1.twiny()
    ax1a.set_yscale('log')
    ax1a.set_xscale('log')
    ax1a.set_xticks([1./18,1./10,1/5.,LS,1])
    ax1a.set_xticklabels(['','10','','$\mathrm{L_s}$','1'],size=18)  
    ax1a.set_xlim(1/18.,2.5)
    ax1a=ax1.twinx()
    ax1a.set_yscale('log')
    ax1a.set_yticks([1./(24*5),1./48,1/24.,1/12])
    ax1a.set_yticklabels(['','','',''])
    ax1a.set_ylim(1/(24*5.),1.)
    if nf==2:
        ax1a.set_yticks([1./(24*5),1./48,1/24.])
        ax1a.set_yticklabels(['5','2','1'],size=18)
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
    if nf==1: ax1.text(1/16,1/(1.7),'a)',color='black' ,size=24)
    if nf==2: ax1.text(1/16,1/(1.7),'b)',color='black' ,size=24)
    
    if nf==1: ax1.text(1/70,1/(2.5*24),'$\mathrm{\omega}$: Frequency [cph]',color='black',size=22,rotation=90)
    if nf==2: ax1.text(1/0.25,1/(2*24),"Period [days]",color='black',size=22,rotation=-90)    
    if nf==2: ax1.text(1/370,1/(15*24),'$\mathrm{\kappa}$: Horizontal wavenumber [cpkm]',color='black',size=22)
            
    if nf==2: 
        ax1.text(1/100,1/(0.018*24),"Wavelength [km]",color='black',size=22)
        cbar_ax = fig.add_axes([0.645,0.07,0.3,0.03])
        aaa=fig.colorbar(cs2, cax=cbar_ax,
                     norm=mpl.colors.Normalize(),
                     ticks=[-100,-50,0,50,100],orientation = 'horizontal')
        aaa.set_ticklabels(['-100','-50','0','50','100'])
        ax1.text(1/0.85,1/(24*30),r'$\mathrm{\kappa \omega }$ Re $[\mathrm{\rho C_p}$ $\mathrm{\widehat{w}\widehat{\theta_a}]}}$ $~\mathrm{(W/m^2)}$',color='black',size=16,horizontalalignment='center')

    # Show the plot
    plt.show()
#=============================================
dirsave='/Users/oceano/MIT-JOB-TESIS/FIGURES/'
plt.savefig(dirsave+'Fig7paper1CoSpectvesion2.png', dpi=600, facecolor='w', edgecolor='w',
        orientation='portrait', papertype=None, format='png',
        transparent=False, bbox_inches=None, pad_inches=0.1,
        frameon=None, metadata=None)









#===========
import numpy
import scipy.special
import scipy.misc
from enthought.mayavi import mlab

r = lambda x,y,z: numpy.sqrt(x**2+y**2+z**2)
theta = lambda x,y,z: numpy.arccos(z/r(x,y,z))
phi = lambda x,y,z: numpy.arctan(y/x)
#phi = lambda x,y,z: numpy.pi+numpy.select(
#	[x>0, x==0, x<0],
#	[
#		numpy.arctan(y/x),
#		.5*numpy.pi*numpy.sign(y),
#		numpy.arctan(y/x)+numpy.pi*numpy.sign(y)]
#)
a0 = 1.
R = lambda r,n,l: (2*r/n/a0)**l * numpy.exp(-r/n/a0) * scipy.special.genlaguerre(n-l-1,2*l+1)(2*r/n/a0)
WF = lambda r,theta,phi,n,l,m: R(r,n,l) * scipy.special.sph_harm(m,l,phi,theta)
absWF = lambda r,theta,phi,n,l,m: abs(WF(r,theta,phi,n,l,m))**2

x,y,z = numpy.ogrid[-24:24:55j,-24:24:55j,-24:24:55j]

mlab.figure()

#mask = numpy.select([theta(x,y,z)>numpy.pi/3.],[numpy.select([abs(phi(x,y,z))<numpy.pi/3.],[numpy.nan],default=1)],default=1)
mask = 1

for n in range(2,3):
	for l in range(1,n):
		for m in range(-l,l+1,1):
			w = absWF(r(x,y,z),theta(x,y,z),phi(x,y,z),n,l,m)
			mlab.contour3d(w*mask,contours=6,transparent=True)

mlab.colorbar()
mlab.outline()
mlab.show()

#  codigo para 
#cmap =plt.cm.get_cmap('seismic')

# extract all colors from the .jet map
#cmaplist = [cmap(i) for i in range(cmap.N)]
# force the first color entry to be grey
#cmaplist[0] = (.5, .5, .5, 1.0)

#cmaplist[66-1:129-1]=cmaplist[1-1:64-1]

#for i in range(0,64,1):    
#    cmaplist[i]=cmaplist[0]
# create the new map
#cmap = mpl.colors.LinearSegmentedColormap.from_list(
#    'Custom cmap', cmaplist, cmap.N)


