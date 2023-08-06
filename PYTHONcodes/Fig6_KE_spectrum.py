#===================Cargar DATOS======AQUI================
# Detalles are in iMac
#
#
#==========================================================
    dirdata='..DATA/'
    name='Ek1TEST.npz'
    name=datos=np.load(dirdata+name)
    Eiso=datos['Eiso'];Eiso=Eiso[:,:-1];
    kiso=datos['kiso'];
    om=datos['om'];om=om[:-1];
    #=======================================================
    print('::::::::: Dimensions :::::::::')
    print('Kiso:',kiso.shape)
    print('omega: ',om.shape)
    print('E: ',Eiso.shape)
    #============================================================
    mpl.rcParams['axes.linewidth'] = 2.5
    mpl.rc('xtick',labelsize=18)
    mpl.rc('ytick',labelsize=18)
    mpl.rc('text',usetex=False)
    mpl.rcParams['xtick.direction'] = 'out'
    mpl.rcParams['ytick.direction'] = 'out'
    mpl.rcParams['font.family']='Times New Roman'
    LS=1/2.3; Ls=2.3
    letra='abcdefghijklmn'
    #==========PLOT  aqui plot plot ==================================================
    # datos to calcluate 
    cmap=plt.cm.get_cmap('CMRmap_r')
    #=============NOrmalizado==============
    scale=10**(0)
    linf=-12; lsup=-4
    bounds2 = np.linspace(linf,lsup,128*2)
    bounds=np.concatenate((bounds2), axis=None)
    normN = mpl.colors.BoundaryNorm(bounds, cmap.N)
    #================KE in vaenumber. Integrate in frequency=====
    name='Ek1TEST.npz'   ; datos=np.load(dirdata+name); Eiso=datos['Eiso'];Eiso=Eiso[:,:-1];
    kE=(Eiso.T).sum(axis=0) 
    fE=(Eiso.T).sum(axis=1)
    #============================================================
    name='Ek4TEST.npz'   ;datos=np.load(dirdata+name);Eiso=datos['Eiso'];Eiso=Eiso[:,:-1];
    kE4=(Eiso.T).sum(axis=0) 
    fE4=(Eiso.T).sum(axis=1)
    #============================================================
    name='Ek6TEST.npz'   ;datos=np.load(dirdata+name);Eiso=datos['Eiso'];Eiso=Eiso[:,:-1];
    kE6=(Eiso.T).sum(axis=0) 
    fE6=(Eiso.T).sum(axis=1)
    #============================================================

    range4x=[1/Ls, 1/1];
    range4y=[range4x[0]**-(2)/63, range4x[1]**-(2)/63]; 
    range3x=[1/Ls, 1/1];
    range3y=[range3x[0]**-(5/3)/1, range3x[1]**-(5/3)/1];
    range2x=[1/8.7, 1/1];
    range2y=[range2x[0]**-(2)/4, range2x[1]**-(2)/4];
    range1x=[1/8.7, 1/Ls];
    range1y=[range1x[0]**-(3)/150, range1x[1]**-(3)/150];



    plt.close('all')
    fig = plt.figure(figsize=(10.,6.5)); plt.show()
    #==============================================================================
    for nf in [1,2,3,4]:#range(1,15):
        print (nf)
        if nf==1:
            ax1 = plt.subplot2grid((28,28),(0,6),rowspan=10,colspan=10)
        if nf==2:
            ax1 = plt.subplot2grid((28,28),(11,0),rowspan=15,colspan=5)  
        if nf==3:
            ax1 = plt.subplot2grid((28,28),(11,6),rowspan=15,colspan=10)
        if nf==4:
            ax1 = plt.subplot2grid((28,28),(11,17),rowspan=15,colspan=10)
    
  
        if nf==1:
            ax1.loglog([1/8.7,1/8.7],[10**-5.1,10**4],'--',linewidth=1.2,color='silver')
            ax1.loglog([1/Ls,1/Ls],[10**-5.1,10**4],'--',linewidth=1.2,color='silver')
            ax1.loglog(kiso,kE6,'-',linewidth=2.5,color='gray',label='Forced(KPP off)')
            ax1.loglog(kiso,kE4,'-',linewidth=2.5,color='black',label='Forced')
            ax1.loglog(kiso,kE,'--',linewidth=2.5,color='black',label='Unforced')
            #ax1.legend(loc='right', bbox_to_anchor=(0.0, 1.3), fontsize=12)
            #ax1.loglog(kiso,(kiso**-3)/150,'--',linewidth=1.0,color='gray')
            #ax1.loglog(kiso,(kiso**-2)/4,'--',linewidth=1.0,color='gray')
            ax1.loglog(range1x,range1y,'--',linewidth=1.0,color='gray')
            ax1.loglog(range2x,range2y,'--',linewidth=1.0,color='gray')
            ax1.loglog(range3x,range3y,'--',linewidth=1.0,color='gray')
            ax1.loglog(range4x,range4y,'--',linewidth=1.0,color='gray')
            ax1.set_xlim(1/35,1/0.4);ax1.set_ylim(10**-3.0,10**2.5)
            ax1.set_xticks([]) ; ax1.set_xticklabels([])
            ax1a=ax1.twiny()
            ax1a.set_yscale('log') ;ax1a.set_xscale('log')
            ax1a.set_xlim(1/35.,2.5)
            ax1a.set_xticks([1./30,1./10,1/5.,LS,1])
            ax1a.set_xticklabels(['30','10','','$\mathrm{L_s}$','1'])
            #ax1.text(1/2.7,2.5*10**6,"Wavelength [km]",color='black',size=18,horizontalalignment='center')
            ax1.text(1/4.0,1.4*10**0.7,"$\kappa^{-2}$",color='gray',size=15,horizontalalignment='center')
            ax1.text(1/4.0,0.32*10**-1,"$\kappa^{-3}$",color='gray',size=15,horizontalalignment='center')
            ax1.text(1/1.3,0.3*10**-2,"$\kappa^{-2}$",color='gray',size=15,horizontalalignment='center')
            ax1.text(1/0.65,3.0*10**1,(letra[nf-1]+')'),color='black',size=21)        
    #        ax1.set_ylim(10**(0),10**(-5))
    #        ax1.invert_xaxis()
            ax1.set_ylabel(r'$\mathrm{KE}(\mathrm{\kappa})[m^2/s^2/cpkm]$',size=12)  
            ax1a.set_xlabel("Wavelength [km]",size=18) 




        if nf==2:
            #ax1.loglog(fE5,om,linewidth=2.0,color='gray')
            ax1.loglog(fE6,om,'-',linewidth=2.0,color='gray',label='Forced (KPP off)')
            ax1.loglog(fE4,om,'-',linewidth=2.0,color='black',label='Forced')
            ax1.loglog(fE,om,'--',linewidth=2.0,color='black',label='Unforced')
            ax1.legend(loc='right',frameon=False, bbox_to_anchor=(0.8,1.31,0.1,0.7), fontsize=14)
            ax1.loglog([10**(5),10**(-8)],[1/24,1/24],'--',linewidth=1.2,color='silver')
            ax1.loglog(om**-3/25500,om,'--',linewidth=1.0,color='gray')
            ax1.loglog(om**-(2)/25,om,'--',linewidth=1.0,color='gray')
            ax1.set_ylim(1/(24*10),1/1)
            ax1.set_xlim(10**(3),10**(-5))
    #        ax1.invert_xaxis()
            ax1.set_ylabel(r'$\mathrm{\omega}$: Frequency [cph]',size=18)       
            ax1.set_xlabel(r'$\mathrm{KE}(\mathrm{\omega})[m^2/s^2/cph]$',size=12) 
            ax1.text(10**-3.5,1.5*10**-1,"$\omega^{-3}$",color='gray',size=13,horizontalalignment='center')
            ax1.text(10**0.90,2.5*10**-0.65,"$\omega^{-2}$",color='gray',size=13,horizontalalignment='center')
            ax1.text(1*10**-3,0.6*10**-2,(letra[nf-1]+')'),color='black',size=21) 
        cmap=plt.cm.get_cmap('CMRmap_r')
    #    cmap=plt.cm.get_cmap('gist_stern_r')
    
        if nf==3:
            name='Ek1TEST.npz'   ; fun='Unforced'
            datos=np.load(dirdata+name)
            Eiso=datos['Eiso'];Eiso=Eiso[:,:-1];
            Em=Eiso.T*kiso[None,...]*om[...,None]; colo='limegreen'        
        if nf==4:
            name='Ek4TEST.npz' ; fun='Forced';  colo='dodgerblue'
            #name='Ek4.npz' ; fun='Heating'; colo='orangered'
            datos=np.load(dirdata+name)
            Eiso=datos['Eiso'];Eiso=Eiso[:,:-1];
            Em=Eiso.T*kiso[None,...]*om[...,None];  
        #=========================================    
        if nf==3 or nf==4:
            ax1.loglog([(1/47)/3,80/47],[(1/180)/3,150/180],'-',linewidth=1.2,color='dimgray')
            #ax1.loglog([(1/4.7),80/4.7],[(1/24)*1.6,150/24*1.6],'-',linewidth=1.2,color='dimgray')
            ax1.text(1/2,1/3.6,r'$\omega$=$c\kappa$',color='dimgray',size=12,rotation=42)
            #cs=plt.pcolormesh(kiso,om,Em, cmap=cmap,shading='gouraud')         
            cs=plt.pcolormesh(kiso,om,np.log(Em), cmap=cmap,shading='gouraud')#shading='gouraud')              
            ax1.set_yscale('log');ax1.set_xscale('log')
    
            ax1.loglog([1/8.7,1/8.7],[1/(24*10.),1.],'--',linewidth=1.2,color='silver')
            ax1.set_ylim(1/(24*10.),1.)
            ax1.set_xlim(1/35.,2.5)
            ax1a=ax1.twiny()
            ax1a.set_yscale('log');ax1a.set_xscale('log')              
            if nf==3: 
                ax1a.set_xticks([1./30,1./10,1/5.,LS,1])
                ax1a.set_xticklabels(['','','','',''])  
            if nf==4: 
                ax1a.set_xticks([1./30,1./10,1/5.,LS,1])
                ax1a.set_xticklabels(['30','10','','$\mathrm{L_s}$','1'])
                ax1.text(1/3.9,1/0.40,"Wavelength [km]",color='black',size=18,horizontalalignment='center')
            ax1a.set_xlim(1/35.,2.5)
            ax1a=ax1.twinx()
            ax1a.set_yscale('log')
            if nf==3:
                ax1a.set_yticks([1./(24*7),1./48,1/24.,1/12])
                ax1a.set_yticklabels(['','','',''])
            ax1a.set_ylim(1/(24*10.),1.)
            if nf==4:
                ax1a.set_yticks([1./(24*7),1./48,1/24.,1/6, 1/1])
                ax1a.set_yticklabels(['7','2','1','6hr','1hr'],size=18)
            ks = np.array([1.e-3,2.5])
            ks = np.array([1.e-3,2.5])
            f  = 1./17.4
            D1= 1/24.0
            ax1.plot(ks,[f,f],'k--',linewidth=1.2,color='black')
            #ax1.text(1.5,f+0.0025,r'$f$',color='black',size=14)
            ax1.plot(ks,[D1,D1],'k--',linewidth=1.2,color='silver')
            ax1.plot([LS,LS],[1/(24*30),1],'k--',linewidth=1.2,color='silver')    
            #if  nf == 3: 
                #ax1.text(1/4.1,1/(4*24),fun,color=colo,size=18)
            #else: 
                #ax1.text(1/4.1,1/(4*24),fun,color=colo,size=18)
            ax1.text(1/30,1/(1.8),(letra[nf-1]+')'),color='black',size=21)        
            if nf==3 or nf==4:
                ax1.set_yticks([]) ; ax1.set_yticklabels([])
            if nf==11 or nf==12 or nf==13 or nf==14: ax1.set_yticks([]) ; ax1.set_yticklabels([])    
        #============================================================#============================================================
            if nf==3: ax1.text(5,1/(45*24),'$\mathrm{\kappa}$: Horizontal wavenumber [cpkm]',color='black',size=18,horizontalalignment='center')
            if nf==4: ax1a.set_ylabel("Period [days]",size=18,rotation=-90, verticalalignment='bottom',horizontalalignment='center')        
        #============================================================#============================================================  
            clim = [-13,-4]
            #clim = [0,2*10**-3]
            plt.clim(clim)  
            cbar_ax = fig.add_axes([0.62,0.77,0.23,0.02])
            aaa=fig.colorbar(cs, cax=cbar_ax,
                         norm=mpl.colors.Normalize(),
                         orientation = 'horizontal')
            if nf==4:
                ax1.text(1/4.1,1/0.045,r'log($\mathrm{\kappa}$ $\mathrm{\omega}$ $\mathrm{\hat{KE}}$($\mathrm{\kappa}$ ,$\mathrm{\omega}$))',color='black',size=16,horizontalalignment='center')  
         #===============================Save grafico !!!!!!!!!!!!================