from astropy.io import fits
import argparse as ap
import numpy as np
from sklearn import svm
from scipy import integrate
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from flare_energy import *
from copy import *

K_B = 1.3806485e-23
C_PLANCK = 6.62607004e-34
C_C = 299792458.0
SIGMA_SB = 5.670374419e-8
NANOMETER = 1.0e-9
LAMBDA_MIN = 454.06
LAMBDA_MAX = 1129.15


def detrended_curve(time, flux):

    """    
    Fits a support vector machine to the curve (time, flux defined)
    with radial kernel functions; returns the trended curve.
    requires the line "from sklearn import svm"
    """

    meant = np.mean(time)
    stdt = np.std(time)
    meanf = np.mean(flux)
    stdf = np.std(flux)

    wt = [[(t - meant) / stdt] for t in time]  # normalization for optimal SVM
    clf = svm.SVR(kernel="rbf", gamma="auto")
    clf.fit(wt, (flux - meanf) / stdf)
    flux_pred = clf.predict(wt)

    return flux_pred * stdf + meanf


def flare_cand_spot(dt_flux, sigma):

    """
	Returns the indices of the flare candidates within a window.
	Ignores nans.
	Goes through jumps, and tags them if they appear significant.
	"""

    flare_times = np.zeros(len(dt_flux))
    len_flare = 0
    df = [dt_flux[i + 1] - dt_flux[i] for i in range(len(dt_flux)-1)]

    # Tag the points where we spot the flare getting over 2 sigma;
    # two consecutive points will suffice:

    for i in range(len(dt_flux)):
        if dt_flux[i] > 3.0 * sigma:
            len_flare += 1
        elif len_flare > 1 and dt_flux[i] > 1.0 * sigma:
            len_flare += 1
        elif len_flare > 1:
            flare_times[i - len_flare : i] = 1.0
            len_flare = 0
        else:
            len_flare = 0

    return flare_times

def flare_ruleout(dt_flux, flare_times):
    
    """
    Rules out not-so-good flares, and also irons out the flare detection pipeline
    """
    
    return np.zeros(len(dt_flux))

def trend_lightcurve(time_raw, pdcflux_raw, windowsize, windowstep):

    """
    the trend of the lightcurve - up to 7 overlapping windows fit the curve with
    svm, take the median of them
    includes NaNs!
    """

    flux_raw = deepcopy(pdcflux_raw)

    n_votes = np.zeros(len(time_raw))
    lightcurve_trended = np.zeros(len(time_raw))
    flare_points = np.zeros(len(time_raw))

    for i in range(0, len(time_raw), windowstep):
        iend = min(len(time_raw), i + windowsize)
        windowflux_raw = deepcopy(flux_raw[i:iend])
        windowtime_raw = time_raw[i:iend]
        has_errs = np.isnan(windowflux_raw)
        flux = windowflux_raw[~has_errs]
        time = windowtime_raw[~has_errs]
        flux_trend = detrended_curve(time, flux)
        windowflux_raw[~has_errs] = flux_trend
        lightcurve_trended[i : i + windowsize] += windowflux_raw
        n_votes[i : i + windowsize] += 1

    lightcurve_trended[np.where(lightcurve_trended == 0)[0]] = np.nan
    lightcurve_trended /= n_votes
    lightcurve_detrended = pdcflux_raw - lightcurve_trended
    errs = np.isnan(lightcurve_detrended)
    sigma = np.std(lightcurve_detrended[~errs])

    for i in range(0, len(time_raw) - windowsize, windowstep):
        
        flare_points[i : i + windowsize] = flare_cand_spot(
            lightcurve_detrended[i : i + windowsize], sigma
        )
    
    return lightcurve_trended, flare_points


def correct_flare(time, pdcflux, fc_start, fc_end, fit_scale, fit_interval):
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
    wt = np.array([[(t-meant)/stdt] for t in time[fc_start-fit_interval-fit_scale : fc_end+fit_interval+fit_scale]])
    
    clf = svm.SVR(kernel="rbf", gamma="auto") 
    clf.fit((ft-meant)/stdt, (ff-meanf)/stdf)
    flux_pred = clf.predict(wt)
    flux_pred = flux_pred*stdf + meanf
    wt = wt*stdt + meant
    plt.plot(wt, pdcflux[fc_start-fit_interval-fit_scale : fc_end+fit_interval+fit_scale], "k.",alpha=0.1)
    plt.plot(wt, flux_pred)
    plt.show()
    

def flare_corrections(time_raw, pdcflux_raw, flaretimes, fit_scale):
    
    """
    Flare times, corrected with a different method from the candidates;
    more reliable when not right next to the edges of the lightcurve
    """
    
    has_errs = np.isnan(pdcflux_raw)
    pdcflux = pdcflux_raw[~has_errs]
    time = time_raw[~has_errs]
    
    in_flare = False
    
    for i in range(len(flaretimes)):
        if flaretimes[i] == 1 and not in_flare:
            in_flare = True
            cand_begin = i
        elif flaretimes[i] == 0 and in_flare:
            cand_end = i
            correct_flare(time, pdcflux, cand_begin, cand_end, fit_scale, int(0.5*fit_scale))
            in_flare = False
            break
        
    

def lcplot(t, f, color=""):
    """
    Simple lightcurve plot function
    """
    plt.plot(t, f, color + ".", alpha=0.1)
    plt.xlabel("Timestamp (days)")
    plt.ylabel(r"SAP flux ($e^-/s$)")


parser = ap.ArgumentParser(description="Trend a lightcurve from a TESS FITS file")
parser.add_argument("inputfile", metavar="if", type=str, help="Input file")
parser.add_argument("period", type=float, help="Period of the star")
parser.add_argument("--plot", help="Plot the lightcurve", action="store_true")
parser.add_argument("--of", type=str, help="Output file for data")

args = parser.parse_args()

hdul = fits.open(args.inputfile)
lcdata = hdul[1].data
pdcflux_raw = lcdata.field(7)
time_raw = lcdata.field(0)
hdul.close()

PERIOD = args.period
WINDOWPERIOD = 0.4
WINDOW_OVERLAP = 9.0

cadence = len(time_raw) / (time_raw[len(time_raw) - 1] - time_raw[0])

windowsize = int(WINDOWPERIOD * PERIOD * cadence)
windowstep = int(windowsize / WINDOW_OVERLAP)

lc_trend, flares = trend_lightcurve(time_raw, pdcflux_raw, windowsize, windowstep)
flare_corrections(time_raw, pdcflux_raw, flares, int(0.05*windowsize))

for i in range(len(time_raw)):
    print("{} {} {} {}".format(time_raw[i], pdcflux_raw[i], lc_trend[i], flares[i]))

if args.plot:
    plotname = args.inputfile[:-4] + "png"
    lcplot(time_raw, pdcflux_raw)
    plt.plot(time_raw[flares == 1.0], pdcflux_raw[flares == 1.0], ".", alpha=0.5)
    plt.plot(time_raw, lc_trend, "r--")
    plt.savefig(plotname)
