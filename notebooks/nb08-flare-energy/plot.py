import numpy as np
from matplotlib import pyplot as plt
data = np.loadtxt("trended_lightcurve_2.out")
 
times = data[:,0]
flarepoints = data[:,3]
trend = data[:,2]
flux = data[:,1]
flare = flarepoints
ft = times[flare == 1.0]
ff = flux[flare == 1.0]
 
plt.plot(times, flux, "k.")
plt.plot(times,trend)
plt.plot(ft,ff,"r.")
plt.show()