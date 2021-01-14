import numpy as np
from matplotlib import pyplot as plt
from scipy.interpolate import interp1d
from scipy.integrate import simps, quad
from astropy.io import fits
import argparse as ap

SIGMA_SB = 5.670374419e-8
LAMBDA_MIN = 454.06
LAMBDA_MAX = 1129.15
NANOMETER = 1.0e-9
ERG = 1.0e-7
R_SUN = 696340000.0
DAY = 86400
T_FLARE = 9000.0

def Response_TESS():
    
    '''
    Returns the TESS response curve as a function of the wavelength
    lambd. It's a spline interpolation 
    of the TESS response function.
    '''
    
    responsefile = "tess-response-function-v1.0.csv"
    rs = np.genfromtxt(responsefile, delimiter=',')
    f = interp1d(rs[:,0], rs[:,1])
    return f

R_tess = Response_TESS()

def B(lambd, T):
    
    '''
    The Planck distribution as a function of wavelength (in nanometers)
    Returns the spectral energy density rather than spectral density.
    '''
    
    l = lambd*NANOMETER
    denom = np.exp(C_PLANCK*C_C/(l*K_B*T)) - 1.
    num = 2. * C_PLANCK * C_C**2. * l**(-5.)
    return num/denom

def L(lambd):
    return B(lambd, T)*R_tess(lambd)

def Lp_STAR(R_star, T_star, R_TESS):

    '''
    Incident luminosity of the star in the wavelength band
    '''
    
    integrand = lambda lam : B(lam, T_star) * R_TESS(lam)
    r = np.arange(LAMBDA_MIN,LAMBDA_MAX)
    integr = simps(integrand(r), r)*NANOMETER
    
    return np.pi * R_star**2 * integr

def Lp_FLARE_DENSITY(T_flare, R_TESS):

    '''
    Incident luminosity density of the flare in the TESS band, divided by its
    area
    '''
    
    integrand = lambda lam : B(lam, T_flare) * R_TESS(lam)
    r = np.arange(LAMBDA_MIN,LAMBDA_MAX)
    integr = simps(integrand(r), r)*NANOMETER
    
    return integr
    
    
def A_flare(rel_amp, LpS, LpF, R_STAR):

    '''
    Returns the estimated effective area of the flare, given the
    relative amplitude from the normalized lightcurve;
    the luminosity of the star, and the luminosity of the flare per unit area
    '''
    
    return rel_amp*LpS/LpF#*R_STAR
    

def L_flare(T_flare, A_f):
    
    '''
    The luminosity of the flare at a given time:
    takes in the area of the flare and its temperature
    '''
    
    return SIGMA_SB * T_flare**4. * A_f
    

def E_flare(amps, times, T_flare, LpS, LpF, R_star):
    
    '''
    The energy of the flare. Takes in the amplitude,
    and pre-calculated values for the luminosity of the 
    star and the luminosity of the flare per unit area
    '''
    
    L_f = lambda amp : L_flare(T_flare, A_flare(amp, LpS, LpF, R_star)) 
    integral = simps(L_f(amps), times)
    return integral



def interflare_time(flare_points, time):

    # Returns a series of interarrival times between flares in days
    
    cadence_time = time[1] - time[0]
    if_times = []
    in_flare = False
    ctr = 0

    for i in range(len(time)):
        if flare_points[i] == 1 and ctr > 1:
            in_flare = True
            if_times.append(cadence_time*ctr)
            ctr = 0
        if flare_points[i] == 1:
            in_flare = True
            ctr = 0
        elif flare_points[i] == 0 and in_flare:
            in_flare = False
            ctr += 1
        else:
            ctr += 1




parse = ap.ArgumentParser(description="Estimate flare energies from a lightcurve with flagged flare times")
parser.add_argument("inputfile", metavar="if", type=str, help="Input file")
parser.add_argument("starinfo, type=str, help = "")
args = parser.parse_args()

# Get information about the star
                    
hdul = fits.open(args.starinfo)
hdr0 = hdul[0].header
T_STAR = hdr0['TEFF']
R_STAR = hdr0['RADIUS']*R_SUN
hdul.close()
                    
LpS = Lp_STAR(R_STAR, T_STAR, R_tess)
LpF = Lp_FLARE_DENSITY(T_FLARE, R_tess)
                    
# Load the lightcurve with flagged flare times

data = np.loadtxt(args.inputfile)
pdcflux = data[:,1]
time = data[:,0]
trend = data[:,2]
flarestamps = data[:,3]
pdcflux_ratio = (pdcflux-trend)/trend

flare_energies = []
flare_times = []
in_flare = False
flare_dur = 0

for i in range(len(pdcflux)):
    if flarestamps[i] == 1:
        in_flare = True
        flare_dur += 1
    elif in_flare:
        amps = pdcflux_ratio[i-flare_dur : i]
        ts = time[i-flare_dur : i]*DAY
        Ef = E_flare(amps, ts, T_FLARE, LpS, LpF, R_STAR)
        flare_energies.append(Ef/ERG)
        flare_times.append(time[i-flare_dur])
        in_flare = False
        flare_dur = 0

# Print the flare energies to stdout

for i in range(len(flare_energies)):
    print("{:.2e} erg, time: {:.1f}".format(flare_energies[i], flare_times[i]))

