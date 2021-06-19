from datetime import datetime, timedelta
import warnings
warnings.filterwarnings("ignore","UserWarning")
warnings.filterwarnings("ignore","RuntimeWarning")
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from matplotlib.lines import Line2D
from netCDF4 import Dataset
import matplotlib.ticker as ticker
import xarray
import scipy.ndimage as ndimage
from scipy.ndimage.filters import gaussian_filter
from scipy.interpolate import make_interp_spline
import glob,os
from scipy import interpolate
import scipy
from statsmodels.nonparametric.smoothers_lowess import lowess
#My script
import ceilometer_help
import peak_finding
import display
from copy import deepcopy
from tqdm import tqdm
from scipy.signal import savgol_filter

read_folder = "/scra6/ft21894/py_script/TAHOPE/obs/PECAN15_SPolCL31/processed/proc2/"
filelist = sorted(glob.glob(read_folder+'*'))
del read_folder

dataCPS = []
for i in range(len(filelist)):
    dataCPS.append(xarray.open_dataset(filelist[i]))

CWTpeak_CPS,backscatterCPS,gpeak_CPS = [],[],[]
for i in range(len(filelist)):
    CWTpeak_CPS.append(dataCPS[i]['CWTpeak_out'])
    backscatterCPS.append(dataCPS[i]['scatterer_out'])
    gpeak_CPS.append(dataCPS[i]['gpeak_out'])

def first_look(CWTarray=None,scatter_array=None,zero_h=50,kdepeak=1,TYPE='morning'):
    tmp = deepcopy(CWTarray)
    tmpc = peak_finding.CWT_transform(scatter_array).clean_peak(peakloc=
                                                                tmp,TYPE='NO_ZERO',zero_loc=zero_h)
    for i in range(0,4):
        NAN = np.isnan(tmpc[:,i])
        tmpc[:,i][NAN] = 0
        del NAN    
    display.plot_CWT_point(CWT_peak=tmpc,plot_peak=[0,1,2],
                           kdepeak=kdepeak,xylim=[0,240*19,0,300],
                           figylim=[0,300],TYPE=TYPE) #CPS08: 0:225*18,0:200
    return tmpc

TIMESTEP = 2
print(filelist[TIMESTEP])

test0802 = first_look(CWTpeak_CPS[TIMESTEP][0:240*19,:],backscatterCPS[TIMESTEP][i,0:300],
                      zero_h=250,kdepeak=1)
display.plot_bkgd(backscatterCPS[TIMESTEP][:,0:300],figylim=[0,300],xylim=None,
                  EXP='PECAN_SPol',cmap='inferno',
                  aspect=6.9,rtitle='McCracken, KS',ltitle='4$^{th}$ June, 2015 (-5)')

display.plot_bkgd(backscatterCPS[TIMESTEP][:,0:300],xylim=[0,240*19,0,300],EXP='PECAN_SPol',cmap='inferno',
                  aspect=6.9,rtitle='McCracken, KS',ltitle='4$^{th}$ June, 2015 (-5)',withline=True,
                  line=[ceilometer_help.smooth(np.ma.masked_equal(test0802[:,0],0),30),
                        ceilometer_help.smooth(np.ma.masked_equal(test0802[:,1],0),30)-9,
                        ceilometer_help.smooth(np.ma.masked_equal(test0802[:,2],0),30)-9],
                  figylim=[0,300],numpeak=3) 
          
def output_KDEh(cleanedpeak=None,peakindex=None,scatter_array=None,intrpsize=None,plt_value=3e-6):
    for num,index in enumerate(peakindex):
        fdisttmp,CWTtmp,CWTtmpintrp = \
            peak_finding.CWT_transform(scatter_array).KDE_manual(xinfo=
                                                                 [0, 240*19, 200j,200],yinfo=[0, 300, 300j,300],
                                                    peakloc=cleanedpeak,peaklocINDEX=index)
        fdisttmpintrp = np.zeros((240*19,300))
        for i in range(300):
            fdisttmpintrp[:,i] = np.interp(intrpsize,
                                           np.linspace(0,240*19-1,200),fdisttmp[i,:])
        B_KDE = lowess(CWTtmpintrp,intrpsize,return_sorted=False,
                              is_sorted=True, frac=0.15, it=0)
        B_KDE3e = lowess(peak_finding.peak_from_KDE(fdisttmpintrp).
                                find_index(values=plt_value),
                                intrpsize,return_sorted=False,is_sorted=True, frac=0.15, it=0)
        B_KDE3e_up = lowess(peak_finding.peak_from_KDE(fdisttmpintrp).
                                   find_index(values=plt_value,TYPE='upper'),
                                   intrpsize,return_sorted=False,is_sorted=True, frac=0.15, it=0)
        del fdisttmp,CWTtmp,CWTtmpintrp,fdisttmpintrp
        return B_KDE,B_KDE3e,B_KDE3e_up
    
test1,test2,test3 = output_KDEh(test0802,[1],backscatterCPS[TIMESTEP][i,0:300],
                                np.linspace(0,240*19-1,240*19))
display.plot_CWT_withbkgd(bkgd=np.asarray(backscatterCPS[TIMESTEP][0:240*19,0:300]),
                          gradientpeak=gpeak_CPS[TIMESTEP][0:240*19,0:300],plot_gpeak=2,xylim=[0,240*19,0,300],
                          CWTpeak=[test1-12,test2-12,test3-12],aspect=7.9,figylim=[0,300],
                          cbarflag=True,ltitle='4$^{th}$ June, 2015 (-5)',
                          rtitle='McCracken, KS',TYPE='morning',scatter_minus=0,EXP='PECAN_SPol')

