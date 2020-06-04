import numpy as np
from sklearn import svm

def detrended_curve(time, flux):

    """    
    Fits a support vector machine to the curve (time, flux defined)
    with radial kernel functions; returns the detrended curve.
    requires the line "from sklearn import svm"
    """
    
    meant = np.mean(time)
    stdt = np.std(time)
    meanf = np.mean(flux)
    stdf = np.std(flux)
    
    wt = [[(t-meant)/stdt] for t in time] # normalization for optimal SVM
    clf = svm.SVR()
    clf.fit(wt,(flux-meanf)/stdf)
    flux_pred = clf.predict(wt)
    
    return flux_pred*stdf+meanf


def find_flare_indices(flux, min_len, sens, start_pos):

    """
    Returns a list of lists of index positions for each flare in a detrended 
    lightcurve segment.
    
    flux is necessarily detrended,
    min_len is the minimum duration of a flare in cadences,
    sens is the minimum (eg 3 sigma)
    """
    
    devlen = 0
    in_flare = False
    flares = []
    ongoing_flare = []

    for i in range(len(detrend_flux)):
    
        if (detrend_flux[i] > 3.*sigma):
            in_flare = True
            ongoing_flare.append(i+start_pos)
            devlen += 1
        elif (in_flare):
            in_flare = False
            if (devlen >= FLARE_MIN_LEN):
                flares.append(ongoing_flare)
            ongoing_flare = []
            devlen = 0
    
    return flares
    

def trend_lightcurve(time_raw,flux_raw,windowsize,windowstep):
    
    """
    the trend of the lightcurve - up to 7 overlapping windows fit the curve with
    svm, take the median of them
    includes NaNs!
    """
    
    n_votes = np.zeros(len(time_raw))
    lightcurve_trended = np.zeros(len(time_raw))
    
    for i in range(0,len(time_raw)-windowsize,windowstep):
        windowflux_raw = flux_raw[i:i+windowsize]
        windowtime_raw = time_raw[i:i+windowsize]
        has_errs = np.isnan(windowflux_raw)
        flux = windowflux_raw[~has_errs]
        time = windowtime_raw[~has_errs]
        flux_trend = detrended_curve(time,flux)
        
        windowflux_raw[~has_errs] = flux_trend
        lightcurve_trended[i:i+windowsize] = windowflux_raw
        n_votes[i:i+windowsize] += 1
     
     for i in range(len(lightcurve_trended)):
        if not np.isnan(lightcurve_trended[i]):
            lightcurve_trended[i] /= float(n_votes[i])
     
     return lightcurve_trended


def flare_overlap(flare1,flaredata):

    #TODO works with an array of flare_cand objects
    
    """
    Check if two flare candidates overlap sufficiently to ignore
    Return -1 if there is no match, return index of matching flare if there
    is
    """
    
    fstart = flare1[0]
    fend = flare1[-1]
    for i in range(len(flaredata)):
    
        f = flaredata[i]
        overlap_start = np.max([fstart,f.i_start])
        overlap_end = np.min([fend,f.i_end])
        total_dur = np.max([fend,f.i_end]) - np.min([fstart,f.i_start])
        overlap_portion = float(overlap_end-overlap_start)/float(total_dur)
        
        if overlap_portion > 0.5:
            return i
    
    return -1

def find_flare_energy(time, flux, flare_pos):

    # TODO
    
    return 0.0


