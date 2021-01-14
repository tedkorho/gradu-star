import argparse as ap
import numpy as np
from sklearn import svm
from matplotlib import pyplot as plt

def needs_correction(time, flux, trend, start, end):
    
    """
    Returns true if the flare is judged to have issues with the fit
    Principle: if the flux appears to be under trend for a time prior/after the flare,
    or the trend looks nothing like a decreasing exponential (Bayesian approach: no trend works
    better than exponential), we need to correct it.
    """
    
    if (start - end == 1): return False
    
    return True 

def correct_flare(time, pdcflux, flares, fc_start, fc_end, fit_scale, fit_interval, sigma):
    
    """
    Corrects a flare
    """
    
    flare_dur = fc_end-fc_start
    fit_flux = []
    fit_time = []
    fit_flux.extend(pdcflux[fc_start-fit_interval-fit_scale:fc_start-fit_interval])
    fit_time.extend(time[fc_start-fit_interval-fit_scale:fc_start-fit_interval])
    fit_flux.extend(pdcflux[fc_end+fit_interval:fc_end+fit_interval+fit_scale])
    fit_time.extend(time[fc_end+fit_interval:fc_end+fit_interval+fit_scale])
    ft = np.array(fit_time)
    ff = np.array(fit_flux)
    
    meanf = np.mean(ff)
    stdf = np.std(ff)
    meant = np.mean(ft)
    stdt = np.std(ft)
    
    ft = np.array([[(t - meant) / stdt] for t in ft])
    ff = np.array((ff-meanf)/stdf)
    
    wt = np.array([[(t-meant)/stdt] for t in time[fc_start-fit_interval-fit_scale : fc_end+fit_interval+fit_scale]])
    
    
    clf = svm.SVR(kernel="rbf", gamma="auto") 
    clf.fit(ft, ff)
    flux_pred = clf.predict(wt) 
    
    flux_pred = flux_pred*stdf + meanf
    wt = wt*stdt + meant
    
    len_flare = 0
    
    for i in range(flare_dur+2*fit_interval):
        flares[i+fc_start-fit_interval] = 0;
        if (pdcflux[i + fc_start - fit_interval] - flux_pred[i] > 3.0 * sigma):
            len_flare += 1
        elif (pdcflux[i + fc_start - fit_interval] - flux_pred[i] > 2.0 * sigma and len_flare > 0):
            len_flare += 1
        elif (len_flare > 2):
            flares[i+fc_start-fit_interval-len_flare : i+fc_start-fit_interval] = 1
            len_flare = 0
        else:
            len_flare = 0
    
    wf = pdcflux[fc_start-fit_interval-fit_scale : fc_end+fit_interval+fit_scale]
    wff = flares[fc_start-fit_interval-fit_scale : fc_end+fit_interval+fit_scale]
    plt.plot(wt, wf, "k.",alpha=0.1)
    plt.plot(wt, flux_pred)
    plt.plot(wt[wff == 1], wf[wff == 1], "r.")
    plt.show()
    
    return flares  
            
    
    

def flare_corrections(time_raw, pdcflux_raw, trend, flaretimes, fit_scale, sigma):
    
    """
    Flare times, corrected with a different method from the candidates;
    more reliable when not right next to the edges of the lightcurve
    """
    
    has_errs = np.isnan(pdcflux_raw)
    flux = pdcflux_raw[~has_errs]
    time = time_raw[~has_errs]
    flare = flaretimes[~has_errs]
    
    in_flare = False
    
    for i in range(len(times)):
        if flare[i] == 1 and not in_flare:
            in_flare = True
            cand_begin = i
        elif flare[i] != 1 and in_flare:
            cand_end = i
            if (needs_correction(time, flux, trend, cand_begin, cand_end)):
                flare = correct_flare(time, flux, flare, cand_begin, cand_end, fit_scale, int(0.5*fit_scale), sigma)
            in_flare = False
        elif times[i] == 0:
            in_flare = False
    
    flaretimes[~has_errs] = flare
    return flaretimes


parser = ap.ArgumentParser()
parser.add_argument("inputfile",type=str)
parser.add_argument("sens",type=float)
parser.add_argument("--of",type=str)
args = parser.parse_args()
data = np.loadtxt(args.inputfile)

times = data[:,0]
flux = data[:,1]
trend = data[:,2]
flares = data[:,3]
has_errs = np.isnan(flares)
sigma = np.std(flux[~has_errs] - trend[~has_errs])
noiselvl = sigma*args.sens

newflares = flare_corrections(times, flux, trend, flares, noiselvl)
for i in range(len(time)):
    print("{:.5f} {:.5f} {:.5f} {:d}".format(times[i], flux[i], trend[i], flares[i]))

plt.plot()