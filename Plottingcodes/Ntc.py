import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import os
import sys

import plotstyle






N = [3, 4, 5, 6, 7, 8, 9, 10, 11]

tc = [0.5, 1.275, 2.775, 4.725, 6.125, 9.025, 11.775, 16.45, 20.25]   #t where zz correlations become negative
dtc = [0.05, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.05, 0.05]


NPBC = [6, 8, 10, 12]

tcPBC = [13.25, 6.55, 9.675, 18.5]
dtcPBC = [0.05, 0.05, 0.075, 0.1]


plt.semilogy(N, tc, '.-', label='OBC')
#plt.plot(NPBC, tcPBC, '.-', label='PBC')
plt.xlabel('Sites')
plt.ylabel(r'$t_c$')
plt.legend()
plt.show()
