from scipy.ndimage.filters import gaussian_filter
import warnings
warnings.filterwarnings("ignore","UserWarning")
warnings.filterwarnings("ignore","RuntimeWarning")
import matplotlib.patheffects as PathEffects
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import matplotlib.ticker as mticker
import matplotlib
import pyart.graph.cm as pyartcm
from matplotlib.axes import Axes
import ceilometer_help

def plot_grd(bkgrd=None,gradient_peak=None,xylim=None,plot_peak=[1,2],figylim=[50,200],TYPE='morning',cmap='inferno'):
    #Preparatory settings
    sns.set_style('white')
    sns.set_context('poster',font_scale=0.6)
    prop = fm.FontProperties(fname='/scra6/ft21894/anaconda3/font/avenir next/AvenirNextLTPro-Regular.otf')

    if (TYPE=='morning'):
        if xylim==None:
            xylim=[0,720*8,0,600]
    elif (TYPE=='night'):
        if xylim==None:
            xylim=[720*18,720*24,0,600]
    print(xylim)
    fig = plt.figure(figsize=(8.55,4.55))
    ax = fig.add_subplot(1,1,1)
        
    color = ['#ff94ab','w','#fffa5e','#ffd166','#ff8fdb']
    
    plt.pcolormesh(np.log10(bkgrd)[:,xylim[2]:xylim[3]].transpose(),cmap=cmap,vmin=-4,vmax=-1.75,zorder=0)
        
    for i in range(np.size(plot_peak)):
        plt.plot(np.linspace(xylim[0],xylim[1]-1,xylim[1]-xylim[0])[:],          
                 ceilometer_help.smooth(gaussian_filter(np.ma.masked_greater(gradient_peak[:,plot_peak[i]],300),0.1)[:],80),
                 lw=2,alpha=0.79,c=color[i])
    ax.set_yticks([0,100,200,300,400,500,600])
    ax.set_yticklabels([0,0.5,1,1.5,2,2.5,3],fontproperties=prop,fontsize=13.5)
    ax.set_xticks([0,720*3,720*6,720*9,720*12,720*15,720*18,720*21,720*24])
    ax.set_xticklabels([0,3,6,9,12,15,18,21,24],fontproperties=prop,fontsize=13.5)    
    ax.tick_params(direction='out',length=6,width=3,colors='k',right=True,left=True,bottom=True)

    plt.xlim(xylim[0],xylim[1])
    plt.ylim(figylim[0],figylim[1])
    plt.show()

def plot_bkgd(bkgrd=None,xylim=None,figylim=[0,150],EXP='CPS08',cmap='inferno',aspect=16.5,
              rtitle=None,ltitle=None,withline=False,line=None,linew=2,numpeak=None):
    cbartick=[-4,-3.5,-3,-2.5,-2,-1.5]
    vmin=-4
    vmax=-1.75
    yticks,yticklbls = [0,50,100,150,200,250,300],[0,0.5,1,1.5,2,2.5,3]
    xticks,xticklbls = [0,240*3,240*6,240*9,240*12,240*15,240*18,240*21,240*24],\
            [0,3,6,9,12,15,18,21,24]
    if (EXP=='CPS08'):# and xylim==None:
        xylim=[0,225*24,0,200]
        yticks,yticklbls = [0,25,50,75,100,125,150,175,200],[0,0.5,1,1.5,2,2.5,3,3.5,4]
        xticks,xticklbls = [0,225*3,225*6,225*9,225*12,225*15,225*18,225*21,225*24],\
            [0,3,6,9,12,15,18,21,24]
    if (EXP=='PECAN_SPol'): #and xylim==None:
        xylim=[0,240*19,0,300]
        vmin,vmax=1,3
        cbartick=[1,1.5,2,2.5,3]
        yticks,yticklbls = [0,50,100,150,200,250,300],[0,0.5,1,1.5,2,2.5,3]
        xticks,xticklbls = [0,240*3,240*6,240*9,240*12,240*15,240*18,240*21,240*24],\
            [0,3,6,9,12,15,18,21,24]
        
        
    sns.set_style('white')
    sns.set_context('poster',font_scale=0.6)
    path = '/scra6/ft21894/anaconda3/font/avenir next/AvenirNextLTPro-Regular.otf'
    prop = fm.FontProperties(fname=path)
    prop = fm.FontProperties(fname='/scra6/ft21894/anaconda3/font/avenir next/AvenirNextLTPro-Regular.otf')
    
    fig = plt.figure(figsize=(10.5,5.5))
    ax = fig.add_subplot(1,1,1)
    plt.imshow(np.log10(bkgrd).transpose(),cmap='inferno',
               vmin=vmin,vmax=vmax,zorder=0,aspect=aspect,origin='lower')
    ax.set_yticks(yticks)
    ax.set_yticklabels(yticklbls,fontproperties=prop,fontsize=13.5)
    ax.set_xticks(xticks)
    ax.set_xticklabels(xticklbls,fontproperties=prop,fontsize=13.5)
    ax.tick_params(direction='out',length=6,width=3,colors='k',right=True,left=True,bottom=True)
    
    
    cbar = plt.colorbar(shrink=0.568,format=mticker.FormatStrFormatter('$10^{%1.1f}$'),
                        ticks=cbartick)
    cbar.ax.tick_params(direction='out',length=6,width=3,colors='k')
    cbar.outline.set_linewidth(2.5)
    cbar.set_label(r'$\beta$ (m$^{-1}$ sr$^{-1}$)',fontproperties=prop,fontsize=14.5)

    ax.autoscale(False)
    if withline==True:
        color = ['r','w','#92d9f7','#ffd166','#ff8fdb']
        numpeak=numpeak
        for i in range(numpeak):
            if color[i]=='w':
                plt.plot(line[i],c=color[i],lw=linew)
            else:
                plt.plot(line[i],c=color[i],lw=linew,path_effects=[PathEffects.withStroke(linewidth=3.45, foreground="k")])
    plt.xlabel('Time (UTC)',fontproperties=prop,fontsize=14)
    plt.ylabel('Height (km)',fontproperties=prop,fontsize=14)
    plt.title(rtitle,loc='right',
              fontproperties=fm.FontProperties(fname='/scra6/ft21894/anaconda3/font/avenir next/AvenirNextLTPro-It.otf'),size=15)
    plt.title(ltitle,loc='left',weight='normal',color='k',    
              fontproperties=fm.FontProperties(fname='/scra6/ft21894/anaconda3/font/avenir next/AvenirNextLTPro-Demi.otf'),size=15)
    plt.xlim(xylim[0],xylim[1])
    plt.ylim(figylim[0],figylim[1])
    plt.tight_layout()
    plt.show()

def plot_CWT_point(CWT_peak=None,xylim=None,plot_peak=[0,1,2],kdepeak=1,figylim=[50,200],TYPE='morning'):
    
    if (TYPE=='morning'):
        if xylim==None:
            xylim=[0,720*8,0,600]
    elif (TYPE=='night'):
        if xylim==None:
            xylim=[720*18,720*24,0,600]

    fig = plt.figure(figsize=(8.55,4.55))
    ax = fig.add_subplot(1,1,1)
    
    color=['m','b','r','g','k']
    for i in range(np.size(plot_peak)):
        plt.plot(np.linspace(xylim[0],xylim[1]-1,xylim[1]-xylim[0])[::10],np.ma.masked_greater(
            ceilometer_help.smooth(CWT_peak[::10,plot_peak[i]],2),300),lw=0,marker='o',\
                 mew=1.85,mfc='w',mec=color[i],ms=4,alpha=0.4)
    sns.kdeplot(np.linspace(xylim[0],xylim[1]-1,xylim[1]-xylim[0])[:], np.ma.masked_invalid(CWT_peak[:,kdepeak]),
                cmap="pyart_Theodore16", shade=True, shade_lowest=False)
    plt.xlim(xylim[0],xylim[1])
    plt.ylim(figylim[0],figylim[1])
    plt.show()

def plot_CWT_withbkgd(bkgd=None,gradientpeak=None,CWTpeak=None,xylim=None,\
                      plot_gpeak=2,figylim=[0,200],cbarflag=True,aspect=6,\
                      cbarlbl=r'$\beta$ (m$^{-1}$ sr$^{-1}$)',ltitle='NTU_CL31',rtitle='4$^{th}$ June 2020',\
                      TYPE='morning',scatter_minus=None,EXP='TAHOPE'):
    #Preparatory settings
    sns.set_style('white')
    sns.set_context('poster',font_scale=0.6)
    prop = fm.FontProperties(fname='/scra6/ft21894/anaconda3/font/avenir next/AvenirNextLTPro-Regular.otf')
        
    if (TYPE=='morning'):
        if xylim==None:
            xylim=[0,720*8,0,600]
    elif (TYPE=='night'):
        if xylim==None:
            xylim=[720*18,720*24,0,600]
    cbartick=[-4,-3.5,-3,-2.5,-2,-1.5]
    if (EXP=='TAHOPE'):
        vmin,vmax=-4,-1.75
        ytick = [0,40,80,120,160,200,240,280,320,360,400,440,480,520,560,600]
        yticklbls = [0,0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2,2.2,2.4,2.6,2.8,3]
        xtick = [0,720*2,720*4,720*6,720*8,720*10,720*12,720*14,720*16,720*18,720*20,720*22,720*24]
        xticklbls = [0,2,4,6,8,10,12,14,16,18,20,22,24]
    elif (EXP=='CPS08'):
        vmin,vmax=-4,-1.75
        ytick = [0,25,50,75,100,125,150,175,200]
        yticklbls = [0,0.5,1,1.5,2,2.5,3,3.5,4]
        xtick = [0,225*3,225*6,225*9,225*12,225*15,225*18,225*21,225*24]
        xticklbls = [0,3,6,9,12,15,18,21,24]
    elif (EXP=='PECAN_SPol'):
        vmin,vmax=1,3
        cbartick=[1,1.5,2,2.5,3]
        ytick,yticklbls = [0,50,100,150,200,250,300],[0,0.5,1,1.5,2,2.5,3]
        xtick,xticklbls = [0,240*3,240*6,240*9,240*12,240*15,240*18,240*21,240*24],\
            [0,3,6,9,12,15,18,21,24]
        
    fig = plt.figure(figsize=(9.5,4.5))
    ax = fig.add_subplot(1,1,1)
    plt.imshow(np.log10(bkgd)[:,xylim[2]:xylim[3]].transpose(),cmap='inferno',
               vmin=vmin,vmax=vmax,zorder=0,aspect=aspect,origin='lower')
    
    ax.set_yticks(ytick)
    ax.set_yticklabels(yticklbls,fontproperties=prop,fontsize=13.5)
    ax.set_xticks(xtick)
    ax.set_xticklabels(xticklbls,fontproperties=prop,fontsize=13.5)
    ax.tick_params(direction='out',length=6,width=3,colors='k',right=True,left=True,bottom=True)
    
    if cbarflag==True:
        cbar = plt.colorbar(shrink=0.78,format=mticker.FormatStrFormatter('$10^{%1.1f}$'),
                            ticks=cbartick)
        cbar.ax.tick_params(direction='out',length=6,width=3,colors='k')
        #cbar.ax.set_yticklabels([0,20,40,60,80,100],fontproperties=prop,fontsize=12.5)
        cbar.outline.set_linewidth(2.5)
        cbar.set_label(r'$\beta$ (m$^{-1}$ sr$^{-1}$)',fontproperties=prop,fontsize=14.5)
    
    ax.autoscale(False)
    plt.plot(np.linspace(xylim[0],xylim[1]-1,gradientpeak[:,plot_gpeak].shape[0])[:],
             ceilometer_help.smooth(gaussian_filter(np.ma.masked_greater(gradientpeak[:,plot_gpeak],300),0.1)[:],40),\
             lw=1.5,ls='--',alpha=0.99,c='#fff94d')
    
    plt.scatter(np.linspace(xylim[0],xylim[1]-1,xylim[1]-xylim[0])[0:int(CWTpeak[0].shape[0]):100], \
                np.asarray(CWTpeak[0])[0:int(CWTpeak[0].shape[0]):100], 
                c = 'r', edgecolor = 'w', s = 40, linewidth='0.75')
    plt.scatter(np.linspace(xylim[0],xylim[1]-1,xylim[1]-xylim[0])[0:int(CWTpeak[0].shape[0]-scatter_minus):100], \
                np.asarray(CWTpeak[1])[0:int(CWTpeak[0].shape[0]-scatter_minus):100], 
                c = 'k', edgecolor = 'w', s = 40, linewidth='0.75',alpha=0.897)
    plt.scatter(np.linspace(xylim[0],xylim[1]-1,xylim[1]-xylim[0])[0:int(CWTpeak[0].shape[0]-scatter_minus):100], \
                np.asarray(CWTpeak[2])[0:int(CWTpeak[0].shape[0]-scatter_minus):100], 
                c = 'k', edgecolor = 'w', s = 40, linewidth='0.75',alpha=0.897)
    plt.xlabel('Time (UTC)',fontproperties=prop,fontsize=14)
    plt.ylabel('Height (km)',fontproperties=prop,fontsize=14)
    plt.title(rtitle,loc='right',
              fontproperties=fm.FontProperties(fname='/scra6/ft21894/anaconda3/font/avenir next/AvenirNextLTPro-It.otf'),size=15)
    plt.title(ltitle,loc='left',weight='normal',color='k',    
              fontproperties=fm.FontProperties(fname=\
                                               '/scra6/ft21894/anaconda3/font/avenir next/AvenirNextLTPro-Demi.otf'),size=15)
    plt.xlim(xylim[0],xylim[1])
    plt.ylim(figylim[0],figylim[1])
    plt.show()

def plot_CWT_withbkgd_smoothmax(bkgd=None,gradientpeak=None,CWTpeak=None,xylim=None,\
                      plot_gpeak=2,figylim=[0,200],cbarflag=True,aspect=6,\
                      cbarlbl=r'$\beta$ (m$^{-1}$ sr$^{-1}$)',ltitle='NTU_CL31',rtitle='4$^{th}$ June 2020',\
                      TYPE='morning',numpeak=2,EXP=None):
    #Preparatory settings
    sns.set_style('white')
    sns.set_context('poster',font_scale=0.6)
    prop = fm.FontProperties(fname='/scra6/ft21894/anaconda3/font/avenir next/AvenirNextLTPro-Regular.otf')
        
    if (TYPE=='morning'):
        if xylim==None:
            xylim=[0,720*8,0,600]
    elif (TYPE=='night'):
        if xylim==None:
            xylim=[720*18,720*24,0,600]
    if (EXP=='TAHOPE'):
        ytick = [0,40,80,120,160,200,240,280,320,360,400,440,480,520,560,600]
        yticklbls = [0,0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2,2.2,2.4,2.6,2.8,3]
        xtick = [0,720*2,720*4,720*6,720*8,720*10,720*12,720*14,720*16,720*18,720*20,720*22,720*24]
        xticklbls = [0,2,4,6,8,10,12,14,16,18,20,22,24]
    elif (EXP=='CPS08'):
        ytick = [0,25,50,75,100,125,150,175,200]
        yticklbls = [0,0.5,1,1.5,2,2.5,3,3.5,4]
        xtick = [0,225*3,225*6,225*9,225*12,225*15,225*18,225*21,225*24]
        xticklbls = [0,3,6,9,12,15,18,21,24]
    elif (EXP=='PECAN_SPol'):
        ytick = [0,50,100,150,200,250,300]
        yticklbls = [0,0.5,1,1.5,2,2.5,3]
        xtick = [0,240*3,240*6,240*9,240*12,240*15,240*18,240*21,240*24]
        xticklbls = [0,3,6,9,12,15,18,21,24]        
            
    fig = plt.figure(figsize=(9.5,4.5))
    ax = fig.add_subplot(1,1,1)
    plt.imshow(np.log10(bkgd)[:,xylim[2]:xylim[3]].transpose(),cmap='inferno',
               vmin=-4,vmax=-1.75,zorder=0,aspect=aspect,origin='lower')
    
    ax.set_yticks(ytick)
    ax.set_yticklabels(yticklbls,fontproperties=prop,fontsize=13.5)
    ax.set_xticks(xtick)
    ax.set_xticklabels(xticklbls,fontproperties=prop,fontsize=13.5)
    ax.tick_params(direction='out',length=6,width=3,colors='k',right=True,left=True,bottom=True)
    
    if cbarflag==True:
        cbar = plt.colorbar(shrink=0.78,format=mticker.FormatStrFormatter('$10^{%1.1f}$'),
                            ticks=[-4,-3.5,-3,-2.5,-2,-1.5])
        cbar.ax.tick_params(direction='out',length=6,width=3,colors='k')
        #cbar.ax.set_yticklabels([0,20,40,60,80,100],fontproperties=prop,fontsize=12.5)
        cbar.outline.set_linewidth(2.5)
        cbar.set_label(r'$\beta$ (m$^{-1}$ sr$^{-1}$)',fontproperties=prop,fontsize=14.5)
    
    ax.autoscale(False)
    plt.plot(np.linspace(xylim[0],xylim[1]-1,gradientpeak[:,plot_gpeak].shape[0])[:],
             ceilometer_help.smooth(gaussian_filter(np.ma.masked_greater(gradientpeak[:,plot_gpeak],300),0.1)[:],40),\
             lw=1.5,ls='--',alpha=0.99,c='#fff94d')
    
    colorlist=['r','b','m','k','y']
    for i in range(numpeak):
        plt.scatter(np.linspace(xylim[0],xylim[1]-1,xylim[1]-xylim[0])[0:int(CWTpeak[i].shape[0]):100],\
                    np.asarray(CWTpeak[i])[0:int(CWTpeak[i].shape[0]):100],\
                    c = colorlist[i], edgecolor = 'w', s = 40, linewidth='0.75')
        
    plt.xlabel('Time (UTC)',fontproperties=prop,fontsize=14)
    plt.ylabel('Height (km)',fontproperties=prop,fontsize=14)
    plt.title(rtitle,loc='right',
              fontproperties=fm.FontProperties(fname='/scra6/ft21894/anaconda3/font/avenir next/AvenirNextLTPro-It.otf'),size=15)
    plt.title(ltitle,loc='left',weight='normal',color='k',    
              fontproperties=fm.FontProperties(fname=\
                                               '/scra6/ft21894/anaconda3/font/avenir next/AvenirNextLTPro-Demi.otf'),size=15)
    plt.xlim(xylim[0],xylim[1])
    plt.ylim(figylim[0],figylim[1])
    plt.show()

