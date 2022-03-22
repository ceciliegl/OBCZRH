import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import os
import sys

from math import comb

import plotstyle

try:
    dir = sys.argv[1]
except:
    print('Give project folder as command line argument.')
    exit(1)


def get_data(filename, variables):
    df = pd.read_csv(filename,\
                    delim_whitespace=True, \
                    engine='python', \
                    names=variables)
    return df
    #using pandas to read the data files

#Loop over all Runs in a folder

runs = [r for r in os.listdir(os.path.expanduser('~/Documents/OBCZRH/Data/' + dir)) if os.path.isdir(os.path.expanduser('~/Documents/OBCZRH/Data/' + dir + '/' + r))]

params = open(os.path.expanduser('~/Documents/OBCZRH/Data/' + dir + '/' + runs[0] + '/parameters.txt'))
Nsites = int(params.readline().split()[-1])
Nh = int(params.readline().split()[-1])
params.close()

t = []
energies = []
zzcorrs = []
xxyycorrs = []

for run in runs:

    #Extract t and add to list
    params = open(os.path.expanduser('~/Documents/OBCZRH/Data/' + dir + '/' + run + '/parameters.txt'))

    for line in params:
        words = line.split()
        if len(words) > 0:
            if words[0] == "Nsites":
                Nsites = int(words[-1])
                break

    params.close()

    params = open(os.path.expanduser('~/Documents/OBCZRH/Data/' + dir + '/' + run + '/parameters.txt'))

    for line in params:
        words = line.split()
        if len(words) > 0:
            if words[0] == "tl":
                tnow = float(words[-1])
                break

    params.close()

    t.append(tnow)

    #Knowing the size of the lattice and the number of holes, you know the number of sectors and their size.

    eigvals = []

    infile = open(os.path.expanduser('~/Documents/OBCZRH/Data/' + dir + '/' + run + '/eigvals.txt'), 'r')

    for line in infile:
        nowvals = [float(i) for i in line.split()[1:]]
        eigvals.append(nowvals)

    energies.append(eigvals)


    infile = open(os.path.expanduser('~/Documents/OBCZRH/Data/' + dir + '/' + run + '/Corr.txt'), 'r')

    beta = []
    time = []

    for line in infile:
        words = line.split()
        beta.append(float(words[0]))
        time.append(float(words[1]))

    beta = np.unique(np.array(beta)) #May be sorting the elements?
    time = np.unique(np.array(time)) #May be sorting the elements?

    Nb = len(beta)
    Nt = len(time)

    infile.close()


    infile = open(os.path.expanduser('~/Documents/OBCZRH/Data/' + dir + '/' + run + '/Corr.txt'), 'r')

    corrzreal = np.zeros((Nsites, Nb, Nt))
    corrzimag = np.zeros((Nsites, Nb, Nt))
    corrpmreal = np.zeros((Nsites, Nb, Nt))
    corrpmimag = np.zeros((Nsites, Nb, Nt))

    for line in infile:
        words = line.split()
        b = float(words[0])
        timenow = float(words[1])

        bind = np.where(beta == b)[0]
        tind = np.where(time == timenow)[0]

        words = words[2:]
        for i in range(Nsites):
            cz = [float(word) for word in words[i][1:-1].split(",")]
            cpm = [float(word) for word in words[Nsites+i][1:-1].split(",")]

            corrzreal[i,bind,tind] = cz[0]
            corrzimag[i,bind,tind] = cz[1]
            corrpmreal[i,bind,tind] = cpm[0]
            corrpmimag[i,bind,tind] = cpm[1]

    zzcorrs.append(corrzreal)
    xxyycorrs.append(corrpmreal)


zzcorrs = np.array(zzcorrs)
xxyycorrs = np.array(xxyycorrs)


t = np.array(t)
indices = np.argsort(t)
t = t[indices]

zzcorrs = zzcorrs[indices,:,:,:]
xxyycorrs = 0.5*xxyycorrs[indices,:,:,:]


#Three site cluster analytical results:
def ZThreeSites(beta, t, J):
    return 4*np.exp(-beta*(J/4.-t)) + 2*np.exp(-beta*(J/4.+2*t)) + np.exp(-beta*(-J/4.-2*t)) + 2*np.exp(-beta*(-J/4.-t)) + 2*np.exp(-beta*(-J/4.+t)) + np.exp(-beta*(-J/4.+2*t))


def S0zS1zThreeSites(beta, t):
    J = -1.
    return (1./12.)*(4.*np.exp(-beta*(J/4.-t)) + 2*np.exp(-beta*(J/4.+2*t)) - np.exp(-beta*(-J/4.-2*t)) - 2*np.exp(-beta*(-J/4.-t)) - 2*np.exp(-beta*(-J/4.+t)) - np.exp(-beta*(-J/4.+2*t)))/ZThreeSites(beta, t, J)


def S0pmS1pmThreeSites(beta, t):
    J = -1.
    return (1./12.)*(-np.exp(-beta*(-J/4.-2*t)) + 2*np.exp(-beta*(-J/4.-t)) - 2*np.exp(-beta*(-J/4.+t)) + np.exp(-beta*(-J/4.+2*t)))/ZThreeSites(beta, t, J)


#print(S0zS1zThreeSites(100,0.4))

#Plot energy evolution as function of t.

"""
#Currently only working for one hole.

plt.figure(1)

Enu0 = np.zeros((len(t), Nsites*1))
Enu1 = np.zeros((len(t), Nsites*(Nsites-1)))
Enu2 = np.zeros((len(t), Nsites*comb(Nsites-1,2)))

for i in range(len(t)):
    Enu0[i,:] = np.array(energies[i][0][:])
    Enu1[i,:] = np.array(energies[i][1][:])
    Enu2[i,:] = np.array(energies[i][2][:])


Enu0 = Enu0[indices,:]
Enu1 = Enu1[indices,:]
Enu2 = Enu2[indices,:]

plt.plot(t, Enu0, '.-', c=plotstyle.folklore)
plt.plot(t, Enu1, '.-', c=plotstyle.lover)
plt.plot(t, Enu2, '.-', c=plotstyle.TS1989)
plt.legend([r"All down", r"One up", r"Two up"])
plt.xlabel(r"$t$")
plt.ylabel(r"Energy")

"""


#Also for each t, plot xy and zz correlations for nearest-neighbour.

plt.figure(2)

for i in range(Nb):
    plt.plot(t, zzcorrs[:,1,i,0], '.-', label=beta[i])
    """if dir == 'OneHole/ThreeSiteCluster':
        plt.plot(t, S0zS1zThreeSites(beta[i], t))"""
plt.xlabel(r"$t$")
plt.ylabel(r"$\langle S_0^z S_1^z \rangle$")
plt.legend()


plt.figure(3)

for i in range(Nb):
    plt.plot(t, xxyycorrs[:,1,i,0], '.-', label=beta[i])
    """if dir == 'OneHole/ThreeSiteCluster':
        plt.plot(t, S0pmS1pmThreeSites(beta[i], t))"""
plt.xlabel(r"$t$")
plt.ylabel(r"$\frac{1}{2}\langle S_0^x S_1^x + S_0^y S_1^y \rangle$")
plt.legend()

plt.show()
