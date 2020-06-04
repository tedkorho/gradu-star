import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import argparse as ap
import pandas as pd
import flaredetect_io as io
import flaredetect_tools as tools

class Flare_cand:
    def __init__(self,confidence,t_start,t_end,energy,significance):
        self.conf = confidence # number 1-n_voters
        self.i_start = t_start # index decided by first detection
        self.i_end = t_end # index decided by first detection
        self.energy = energy # an array of votes
        self.sig = significance # sigmas of peak


def main():
    
    args = io.get_args()
    
    hdul = fits.open(args.inputfile)
    lcdata = hdul[1].data
    pdcflux_raw = lcdata.field(7)
    time_raw = lcdata.field(0)
    
    PERIOD = args.p
    SENS = args.sens
    MIN_LEN = args.minlen
    WINDOWPERIOD = 0.4
    WINDOW_OVERLAP = 7
    
    cadence = len(time_raw)/(time_raw[len(time_raw)-1] - time_raw[0])
    windowsize = int(WINDOWPERIOD*PERIOD*cadence)
    windowstep = int(windowsize/WINDOW_OVERLAP)
    
    # Initialize data structure for flare positions, confidence (1-7)
    # and find sigmas
    # -> has to be extendable with energy estimates
    # pandas?
    
    flaredata = [] # AN ARRAY OF FLARE_CAND OBJECTS
    
    # Loop proper:
    
    for i in range(0,len(time_raw)-windowsize,windowstep):
        
        # Remove NaNs, skip window if there are too many
    
        windowflux_raw = pdcflux_raw[i:i+windowsize]
        windowtime_raw = time_raw[i:i+windowsize]
        has_errs = np.isnan(windowflux_raw)
        flux = windowflux_raw[~has_errs]
        time = windowtime_raw[~has_errs]
        
        if float(len(flux))/float(len(windowflux_raw)) < 0.3:
            continue
    
        # Find flares,  analyze energy estimates
        
        flux_detrend = tools.detrended_curve(time,flux)
        flare_pos = find_flare_indices(flux_detrend, 
                        MIN_LEN, SENS*np.std(flux_detrend), i)
        energies = []
        for j in range(len(flare_pos)):
            energies.append(find_flare_energy(time,flux,flare_pos[j])) # TODO
        
        # Consolidate overlapping flares, vote on detection and energy
        
        for j in range(len(flare_pos)):
            ol = tools.flare_overlap(flare_pos[j],flaredata):
            if (ol != -1):
                flaredata[ol].conf += 1
                flaredata[ol].energy.append(energies[j])
                continue
            else
                flaredata.append(Flare_cand(1,
                                    flare_pos[j][0],
                                    flare_pos[j][-1],
                                    energies[j],
                                    0))
            
            
            
    
    # Write outputs, create plots if demanded
    
    io.write_flares(args.o,flaredata) #TODO


if __name__ == "__main__":
    main()

