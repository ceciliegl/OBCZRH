import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import os
import sys
import scipy

import plotstyle


def tcfit(x, a, b, alpha):
    return a + b*x**alpha

N = [3, 4, 5, 6, 7, 8, 9, 10, 11, 12]

tc = [0.5, 1.275, 2.775, 4.725, 6.125, 9.025, 11.775, 16.45, 20.25, 27.20]   #t where zz correlations become negative
dtc = [0.05, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.05, 0.05, 0.05]


popt, pcov = scipy.optimize.curve_fit(tcfit, N, tc)

print(popt)

NPBC = [6, 8, 10, 12]

tcPBC = [13.25, 6.55, 9.675, 18.5]
dtcPBC = [0.05, 0.05, 0.075, 0.1]


NPBCy = [3, 4, 5, 6, 7, 8, 9, 10, 11, 12]  #Only periodic in y, open in x
tcPBCy = [0.75, 0.925, 3.025, 4.275, 6.975, 9.075, 13.475, 16.725, 23.25, 28.0]
dtcPBCy = [0.05, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.05, 0.05]


plt.plot(N, tc, '.-', label='OBC')
plt.plot(N, tcfit(N, *popt))
plt.plot(NPBC, tcPBC, '.-', label='PBC')
plt.plot(NPBCy, tcPBCy, '.-', label='PBCy')
plt.xlabel('Sites')
plt.ylabel(r'$t_c$')
plt.legend()
plt.show()
