import numpy as np
from sklearn.neighbors import KernelDensity
from scipy import interpolate

def window(size):
	    return np.ones(size)/float(size)

def cbhindex(cbharray=None):
	height = np.linspace(0,1449,1500)*5/1000
	cbh_index = np.zeros_like(np.asarray(cbharray))
	
	for i in range(cbh_index.shape[0]):
		if np.isnan(np.asarray(cbharray)[i])==False:	
			cbh_index[i] = np.abs(height-np.asarray(cbharray)[i]/1000).argmin()
			
		else:
			cbh_index[i] = 0
	return cbh_index

def myhaar(a=40,b=None,inputarray=None):
    """
    Define Haar Wavelet here.
    a: Wavelet scale
    b: Location of Wavelet Centre
    """
    haar = np.zeros_like(inputarray)
    
    if (b<a/2):
        haar[0:b] = 1
        haar[b:int(b+a/2)] = -1
    if ((inputarray.shape[0]-a/2)<a/2):
        haar[int(b-a/2):b] = 1
        haar[b:inputarray.shape[0]] = -1
    else:
        haar[int(b-a/2):b] = 1
        haar[b:int(b+a/2)] = -1
    return haar

def mymxh(a=40,b=None,inputarray=None):
    """
    Define Mexican Hat Wavelet here.
    a: Wavelet scale
    b: Location of Wavelet Centre
    """
    mxh = np.zeros_like(inputarray)
    X = np.linspace(0,mxh.shape[0]-1,mxh.shape[0])
    a=a/2.5
    mxh = 10.5*(np.exp(-((X-b)/a)**2/2)-np.exp(-((X-b)/a)**2/2)*((X-b)/a)**2)/(4*np.pi)
    return mxh

def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth

def kde2D(x, y, bandwidth, xbins=100j, ybins=100j, **kwargs): 
    """Build 2D kernel density estimate (KDE)."""

    # create grid of sample locations (default: 100x100)
    xx, yy = np.mgrid[x.min():x.max():xbins, 
                      y.min():y.max():ybins]

    xy_sample = np.vstack([yy.ravel(), xx.ravel()]).T
    xy_train  = np.vstack([y, x]).T

    kde_skl = KernelDensity(bandwidth=bandwidth, **kwargs)
    kde_skl.fit(xy_train)

    # score_samples() returns the log-likelihood of the samples
    z = np.exp(kde_skl.score_samples(xy_sample))
    return xx, yy, np.reshape(z, xx.shape)

def autocovariance(Xi, N, k, Xs):
    autoCov = 0
    for i in np.arange(0, N-k):
        autoCov += ((Xi[i+k])-Xs)*(Xi[i]-Xs)
    return (1/(N-1))*autoCov

def autocorrelation():
    return autocovariance(Xi, N, k, Xs) / autocovariance(Xi, N, 0, Xs)