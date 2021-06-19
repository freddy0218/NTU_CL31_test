import numpy as np
from scipy.ndimage.filters import gaussian_filter
import ceilometer_help
import display

class gradient_method(object):
    
    def __init__(self, scatterer=None,yaxis=None,vertical_resolution=5):
        self.scatterer=scatterer
        self.yaxis=yaxis
        self.vert_res=vertical_resolution
    
    def gradient_peak(self):
        from scipy.signal import find_peaks
        input_peaks = np.zeros((self.scatterer.shape[0],8))
        
        for i in range(self.scatterer.shape[0]):
            TEMP = np.gradient(gaussian_filter(np.log10(self.scatterer[i,:]),0.1),self.yaxis*self.vert_res)
            TEMP[TEMP>-0.005] = 0 #Filter positive and very small negative oscillations
            peaks, _ = find_peaks(-TEMP)
            if (peaks.shape[0]<8) and (peaks.shape[0]>0):
                input_peaks[i,0:peaks.shape[0]] = peaks
                peaks[peaks.shape[0]:8] = peaks[-1]+10
            elif (peaks.shape[0]==0):
                peaks[0:8] = np.nan
            else:    
                input_peaks[i,:] = peaks[0:8]
        return input_peaks

class CWT_transform(object):
    
    def __init__(self, scatterer=None,vertical_resolution=5):
        self.scatterer=scatterer
        self.vert_res=vertical_resolution
    
    def CWT_peaks(self,a=None,TYPE='haar',peak_distance=30):
        from scipy.ndimage.filters import gaussian_filter
        from scipy.signal import find_peaks
        
        #Gaussian filter
        A = gaussian_filter(np.log10(self.scatterer),0.75)
        #Integrate Wavelet Covariance
        if TYPE=='haar':
            haar_CWT = np.zeros_like(A)
            A[np.isnan(A)] = 0
            
            for i in range(A.shape[0]):
                haar_CWT[i] = np.trapz(A[i:]*ceilometer_help.myhaar(a=a,b=i,
                                                                    inputarray=self.scatterer)[i:],
                                       x=np.linspace(0,A.shape[0]-1,A.shape[0])[i:]*self.vert_res)
        #Find peak locations
        peaksCWT, _ = find_peaks(gaussian_filter(haar_CWT,1.08),distance=peak_distance)
        return haar_CWT,peaksCWT
    
        if TYPE=='mxh':
            mxh_CWT = np.zeros_like(A)
            A[np.isnan(A)] = 0
            
            for i in range(A.shape[0]):
                mxh_CWT[i] = np.trapz(A[i:]*ceilometer_help.mymxh(a=a,b=i,inputarray=self.scatterer)[i:],
                                      x=np.linspace(0,A.shape[0]-1,A.shape[0])[i:]*self.vert_res)
        peaksCWT, _ = find_peaks(gaussian_filter(mxh_CWT,1.08),distance=peak_distance)    
        return mxh_CWT,peaksCWT        
    
    def clean_peak(self,peakloc=None,TYPE=None,zero_loc=None):
        
        CWT0_5copy = np.copy(peakloc)
        
        # Sudden jump adjustment
        if TYPE=='SUDDEN':
            for i in range(peakloc.shape[0]-1):
                if (CWT0_5copy[i+1,0]-CWT0_5copy[i,0]>20): #Sudden Jump
                    CWT0_5copy[i+1,3] = CWT0_5copy[i+1,2]
                    CWT0_5copy[i+1,2] = CWT0_5copy[i+1,1]
                    CWT0_5copy[i+1,1] = CWT0_5copy[i+1,0]
                    CWT0_5copy[i+1,0] = CWT0_5copy[i,0]
        elif TYPE=='NO_ZERO':
            for i in range(peakloc.shape[0]-1):
                if (CWT0_5copy[i,0]>zero_loc):
                    CWT0_5copy[i,3]  = CWT0_5copy[i,2]
                    CWT0_5copy[i,2]  = CWT0_5copy[i,1]
                    CWT0_5copy[i,1]  = CWT0_5copy[i,0]
                    CWT0_5copy[i,0]  = np.nan#np.mean(np.ma.masked_greater(CWT0_5copy,zero_loc))
        return CWT0_5copy
    
    def KDE_manual(self,xinfo=None,yinfo=None,peakloc=None,peaklocINDEX=None):
        import scipy.stats as st
        xmin, xmax = xinfo[0],xinfo[1]#0, 720*8
        ymin, ymax = yinfo[0],yinfo[1]#0, 300
        
        # Peform the kernel density estimate
        yy,xx = np.mgrid[ymin:ymax-1:yinfo[2], xmin:xmax-1:xinfo[2]] #300j,200j
        positions = np.vstack([xx.ravel(), yy.ravel()])
        values = np.vstack([np.linspace(0,peakloc.shape[0]-1,peakloc.shape[0]),
                            np.ma.masked_greater(peakloc[:,peaklocINDEX],300)])
        kernel = st.gaussian_kde(values)
        f = np.reshape(kernel(positions).T, xx.shape)
        
        #KDE maximum location
        argmax_f = []
        for i in range(xinfo[3]):
            argmax_f.append(f[:,i].argmax())
            
        return f,np.asarray(argmax_f),np.interp(np.linspace(0,peakloc.shape[0]-1,peakloc.shape[0]),
                                                np.linspace(0,peakloc.shape[0]-1,xinfo[3]),np.asarray(argmax_f))
    

class peak_from_KDE(object):
    
    def __init__(self,KDE_obj=None):
        self.KDE_obj = KDE_obj
    
    def find_percentileindex(self,percentile=95,TYPE='lower'):
        if TYPE=='lower':
            output_index = []
            for i in range(self.KDE_obj.shape[0]):
                temp = np.percentile(self.KDE_obj[i,:],percentile)
                output_index.append(np.abs(self.KDE_obj[i,0:(self.KDE_obj[i,:]).argmax()]-temp).argmin())
        return output_index
    
    def find_index(self,values=None,TYPE='lower'):
        if TYPE=='lower':
            output_index = []
            for i in range(self.KDE_obj.shape[0]):
                output_index.append(np.abs(self.KDE_obj[i,0:(self.KDE_obj[i,:]).argmax()]-values).argmin())
        elif TYPE=='upper':
            output_index=[]
            for i in range(self.KDE_obj.shape[0]):
                output_index.append(np.abs(self.KDE_obj[i,(self.KDE_obj[i,:]).argmax():]-\
                                           values).argmin()+self.KDE_obj[i,:].argmax())
        return output_index
