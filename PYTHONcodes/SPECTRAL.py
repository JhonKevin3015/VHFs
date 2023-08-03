#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 26 15:43:55 2022

@author: oceano
"""
#A=wpc
#B=tpc

import numpy as np
from scipy import signal

def cospec_ab_function(A,B,d1,d2,d3,windowsclass):
    # we try to use the co-spectrum
    # condition: use Hanning in x direction with periodic boundary!
    #A=wpc
    l1,l2,l3 = A.shape
    N = A.size
    print('before fftn')
    print(l1,l2,l3)
    df1 = 1./(l1*d1)
    df2 = 1./(l2*d2)
    df3 = 1./(l3*d3)
    f1Ny = 1./(2*d1)
    f2Ny = 1./(2*d2)
    f3Ny = 1./(2*d3)
    f1 = np.arange(-f1Ny,f1Ny,df1)
    f2 = np.arange(-f2Ny,f2Ny,df2)
    f3 = np.arange(0,l3/2+1)*df3

    # spectral window
    # first, the spatial window
    if windowsclass=='flat':
        wx = np.matrix(signal.windows.flattop(l1,sym=False)*0+1)  
        wy = np.matrix(signal.windows.flattop(l2,sym=False))
        window_s = np.repeat(np.array(wx.T*wy),l3).reshape(l1,l2,l3)
        # now, the time window
        wt = signal.windows.flattop(l3,sym=False)
        window_t = np.repeat(wt,l1*l2).reshape(l3,l2,l1).T

    if windowsclass=='hann': 
        wx = np.matrix(np.hanning(l1)*0+1)    
        wy = np.matrix(np.hanning(l2))
        window_s = np.repeat(np.array(wx.T*wy),l3).reshape(l1,l2,l3)
        # now, the time window
        wt = np.hanning(l3)
        window_t = np.repeat(wt,l1*l2).reshape(l3,l2,l1).T

    if windowsclass=='kais':
        beta=3.86         
        wx = np.matrix(np.kaiser(l1,beta)*0+1)  
        wy = np.matrix(np.kaiser(l2,beta))
        window_s = np.repeat(np.array(wx.T*wy),l3).reshape(l1,l2,l3)
        # now, the time window
        wt = np.kaiser(l3,beta)
        window_t = np.repeat(wt,l1*l2).reshape(l3,l2,l1).T    

    Ahat = np.fft.rfftn(window_s*window_t*A)
    Bhat = np.fft.rfftn(window_s*window_t*B)

    Aabs = (np.abs(Ahat)**2)
    ## Parseval's for A
    e = np.sum(A**2)
    spec_e = np.sum(Aabs)/N
    ratio_eA = e/spec_e

    Babs = (np.abs(Bhat)**2)
    ## Parseval's for B
    e = np.sum(B**2)
    spec_e = np.sum(Babs)/N
    ratio_eB = e/spec_e
    print(" Correction Factor to A and  B matrix: ",ratio_eA*ratio_eB)  
    # compute co-spectrum with Parseval correction factor
    #cospec =(Bhat*Ahat.conjugate()).real / ((l1*l2*l3)**2) / (df1*df2*df3) * np.sqrt(ratio_eA)*np.sqrt(ratio_eB)
    cospec =(Bhat*Ahat.conjugate()).real/ ((l1*l2*l3)**2) / (df1*df2*df3) 
    cospec_dens = np.fft.fftshift(cospec.copy(),axes=(0,1))
    cospec_rms = np.sqrt((cospec_dens[1:,1:,1:-2].sum()*df1*df2*df3))
    return np.fft.fftshift(ratio_eA*ratio_eB*cospec,axes=(0,1)),f1,f2,f3,df1,df2,df3

#============================================================================================
#ver=window_s*A
#plt.figure ()
#plt.pcolor(LX,LX,ver[:,:,300]*24*60*60,cmap=plt.cm.get_cmap('seismic'),shading='flat')
#clim = [-100,100];plt.clim(clim);plt.colorbar(cs,label=r'');ax1.set_title ('V. velocity`')
#ax1.set_aspect('equal', 'box')

#=======This is used to KE and other variable================================================
def spec_est3(A,d1,d2,d3,varname,beta=3.86):
    l1,l2,l3 = A.shape
    N = A.size
    df1 = 1./(l1*d1)
    df2 = 1./(l2*d2)
    df3 = 1./(l3*d3)
    #f1Ny = 1./(2*d1)
    #f2Ny = 1./(2*d2)
    #f3Ny = 1./(2*d3)
    #f1 = np.arange(-f1Ny,f1Ny,df1)
    #f2 = np.arange(-f2Ny,f2Ny,df2)
    #f3 = np.arange(0,l3/2+1)*df3
    f1 = np.fft.fftshift(np.fft.fftfreq(l1,d1))
    f2 = np.fft.fftshift(np.fft.fftfreq(l2,d2))
    f3 = np.fft.rfftfreq(l3,d3)

    # spectral window -- Kaiser window
    # https://dsp.stackexchange.com/questions/40598/why-would-one-use-a-hann-or-bartlett-window
    # https://docs.scipy.org/doc/numpy/reference/generated/numpy.kaiser.html
    # first, the spatial window
    #beta = 4.86
    wx = np.matrix(np.kaiser(l1,beta)*0+1) # le cambiamos por ser un cannal
    wy = np.matrix(np.kaiser(l2,beta))
    window_s = np.repeat(np.array(wx.T*wy),l3).reshape(l1,l2,l3)
    # now, the time window
    wt = np.kaiser(l3,beta)
    window_t = np.repeat(wt,l1*l2).reshape(l3,l2,l1).T
    Ahat = np.fft.rfftn(window_s*window_t*A)
    Aabs = (np.abs(Ahat)**2) #/ (df1*df2*df3) / ((l1*l2*l3)**2)

    ## Parseval's
    # Energy in original domain
    e = np.sum(A**2)#*d1*d2*d3
    print(varname+" Energy (x,y,t)",e)
    #Energy in FFT domain
    spec_e = np.sum(Aabs)/N #*df1*df2*df3
    print(varname+" Energy (k,l,f)",spec_e)
    # Energy ratio
    ratio_e = e/spec_e
    # Correting
    Aabs = Aabs*ratio_e
    print(varname+" Corrected energy (k,l,f)",np.sum(Aabs)/N)

    ## PSD
    print("Dividing by N (total # of elems) to get (approximate) PSD")
    Aabs = Aabs/N

    return np.fft.fftshift(Aabs,axes=(0,1)),f1,f2,f3
#============================================================================================
#:::::::::::::::::::::::::::::::::::::::::::::::
from numpy import  pi
def calc_ispec(k,l,E):
    """ calculates isotropic spectrum from 2D spectrum """

    dk,dl = k[1,]-k[0],l[1,]-l[0]
    l,k = np.meshgrid(l,k)
    wv = np.sqrt(k**2 + l**2)

    if k.max()>l.max():
        kmax = l.max()
    else:
        kmax = k.max()

    # create radial wavenumber
    dkr = np.sqrt(dk**2 + dl**2)
    kr =  np.arange(dkr/2.,kmax+dkr,dkr)
    ispec = np.zeros(kr.size)
    #print(ispec.shape)
    #print(kr.shape)
    for i in range(kr.size):
        fkr =  (wv>=kr[i]-dkr/2) & (wv<=kr[i]+dkr/2)
    #    print(fkr.shape)
        dth = pi / (fkr.sum()-1)
        ispec[i] = E[fkr].sum() * kr[i] * dth

    return kr, ispec
#:::::::::::::::::::::::::::::::::::::::::
    #:::::::::::::::::::::::::::::::::::::::::::
# Gradients of U
# definition: uv vel-comp; d horizontal interval
def gradu(u,v,d):
    ux,uy = np.gradient(u,d,d)
    vx,vy = np.gradient(v,d,d)
    vort = vx - uy
    div = ux + vy
    strain = ((ux-vy)**2 + (vx+uy)**2)**.5
    #ss=vx+uy;
    #sn=ux-vy;
    ow=strain**2-vort**2;        
    return vort,div,ow
#:::::::::::::::::::::::::::::::::::::::::::
import numpy as np
def timtim(ntpos,NT,dt):
    nt=NT[960-1+ntpos]
    tim=(nt+1+960)*dt/(60*60*24)
    dia=np.fix(tim)
    hora=np.fix(np.fix((tim-np.fix(tim))*100)*24/100)
    titex=[str(np.intp(dia))+'d:'+str(np.intp(hora))+'hrs']       
    return titex
#:::::::::::::::::::::::::::::::::::::::::::