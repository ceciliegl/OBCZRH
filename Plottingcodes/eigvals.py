#coding: UTF-8

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import os
import sys

import plotstyle

try:
    run_numbers = sys.argv[1:]
except:
    print("Give run number as command line argument.")
    exit(1)


def get_data(filename, variables):
    df = pd.read_csv(filename,\
                    delim_whitespace=True, \
                    engine='python', \
                    names=variables)
    return df
    #using pandas to read the data files

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


    #data = get_data("Run" + run_number + "/eigvals.txt", ["", "OP1", "OP2", "OP1norm", "OP2norm"])

    eigvals = []

    infile = open("Run" + run_number + "/eigvals.txt", "r")

    for line in infile:
        eigvals += line.split()[1:]

    eigvals = [float(i) for i in eigvals]

    print(np.sort(eigvals)[:100])


    plt.figure()
    plt.hist(eigvals, bins=100)
    plt.xlabel(r"E", fontsize=14)
    plt.ylabel(r"Count", fontsize=14)
    plt.title("Run"+ run_number)
    #plt.title("One hole, Ising FM, t = 100")
plt.show()
