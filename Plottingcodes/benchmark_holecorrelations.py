import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import os
import sys

from math import comb

import plotstyle

try:
    run_numbers = sys.argv[1:]
except:
    print("Give run number as command line argument.")
    exit(1)


for run_number in run_numbers:
    if float(run_number) < 10:
        run_number = "00" + run_number
    elif float(run_number) < 100:
        run_number = "0" + run_number
    elif float(run_number) < 1000:
        run_number = run_number
    else:
        print("Run number too big")
        exit(1)

    Zfile = open("Run" + run_number + "/Partition.txt", "r")

    Zbeta = []
    Z = []

    for line in Zfile:
        words = line.split()
        Zbeta.append(float(words[0]))
        Z.append(float(words[1]))

    Zbeta = np.array(Zbeta)
    Z = np.array(Z)

    Zfile.close()

    paramsfile = open("Run" + run_number + "/parameters.txt", "r")

    for line in paramsfile:
        words = line.split()
        if len(words) > 0:
            if words[0] == "Nsites":
                Nsites = int(words[-1])
                break

    print("Nsites = %d" % Nsites)

    paramsfile.close()

    infile = open("Run" + run_number + "/HoleCorr.txt", "r")

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


    infileZERO = open("Run" + run_number + "/HoleCorr.txt", "r")
    infileNN = open("Run" + run_number + "/HoleCorrNN.txt", "r")
    infileMID = open("Run" + run_number + "/HoleCorrMID.txt", "r")

    corrrealZERO = np.zeros((Nsites, Nb, Nt))
    corrimagZERO = np.zeros((Nsites, Nb, Nt))

    corrrealNN = np.zeros((Nsites, Nb, Nt))
    corrimagNN = np.zeros((Nsites, Nb, Nt))

    corrrealMID = np.zeros((Nsites, Nb, Nt))
    corrimagMID = np.zeros((Nsites, Nb, Nt))

    for line in infileZERO:
        words = line.split()
        b = float(words[0])
        t = float(words[1])

        bind = np.where(beta == b)[0]
        tind = np.where(time == t)[0]

        words = words[2:]
        for i in range(Nsites):
            c = [float(word) for word in words[i][1:-1].split(",")]

            corrrealZERO[i,bind,tind] = c[0]
            corrimagZERO[i,bind,tind] = c[1]

    for line in infileNN:
        words = line.split()
        b = float(words[0])
        t = float(words[1])

        bind = np.where(beta == b)[0]
        tind = np.where(time == t)[0]

        words = words[2:]
        for i in range(Nsites):
            c = [float(word) for word in words[i][1:-1].split(",")]

            corrrealNN[i,bind,tind] = c[0]
            corrimagNN[i,bind,tind] = c[1]

    for line in infileMID:
        words = line.split()
        b = float(words[0])
        t = float(words[1])

        bind = np.where(beta == b)[0]
        tind = np.where(time == t)[0]

        words = words[2:]
        for i in range(Nsites):
            c = [float(word) for word in words[i][1:-1].split(",")]

            corrrealMID[i,bind,tind] = c[0]
            corrimagMID[i,bind,tind] = c[1]



    indices = np.argsort((np.arange(Nsites)+np.ones(Nsites)*Nsites/2)%Nsites)

    plt.figure()
    for b in range(len(beta)):
        plt.plot(np.arange(Nsites), corrrealZERO[:,b,0]-corrrealMID[:,b,0][indices], 'o', label=r"$\beta = %.2f$" % beta[b])
    plt.xlabel(r"site $j$")
    plt.ylabel(r"$\langle N^h_{0}N^h_{j} \rangle$")
    plt.title("ZERO VS MID Run"+ run_number)
    #plt.title("One hole, Ising FM, t = 10")
    plt.legend()


    plt.figure()
    for b in range(len(beta)):
        plt.plot(np.arange(Nsites), corrrealNN[:,b,0]-corrrealNN[0,b,0]*np.ones(Nsites), 'o', label=r"$\beta = %.2f$" % beta[b])
    plt.xlabel(r"site $j$")
    plt.ylabel(r"$\langle N^h_{0}N^h_{j} \rangle$")
    plt.title("NN Run"+ run_number)
    #plt.title("One hole, Ising FM, t = 10")
    plt.legend()



plt.show()
