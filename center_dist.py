# LW 6/11/15
# This script loads transformed IPD as a numpy array, calculates the mean IPDt, and 
# subtracts that mean from each IPD and writes it out to a new file.
# This is to center the IPDt distribution at zero.

#Usage: python center_dist.py

from scipy import stats
import numpy as np
import statistics


f1 = np.loadtxt('24hr_IPDt_less1000.txt', dtype=float, skiprows=1)
mean = np.mean(f1, axis=0)
newIPD = f1[0:] - mean
np.savetxt('24hr_centeredIPD_less1000.txt', newIPD, delimiter=' ', header='newIPDt')
