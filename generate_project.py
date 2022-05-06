import os
import numpy as np

<<<<<<< HEAD
mainproject = "Pyrochlore"  #Set to zero if only one project.
project = "2holes3tetrahedraADDITIONAL"
description = "Testing."
jobname = "myjob"
time = "5:00:00"
runmin = 17
runmax = 21
=======
mainproject = "TwoHolesMIDCORRBettervaluesoft"  #Set to zero if only one project.
project = "N7"
description = "Testing."
jobname = "myjob"
time = "5:00:00"
runmin = 0
runmax = 19
>>>>>>> ef3c47272fadb627936e0457744360c245e4d47c
runsame = 0
nruns = (runmax-runmin) + 1
NICE = 11

#LIBS#
#BOOST = 0       #Higher precision in Eigen-calculations. Time-consuming. Not implemented for now.

#LATTICE#
Nsites = 10*np.ones(nruns, int)
nruns  = len(Nsites)
runmax = runmin + (nruns-1)

Nh = 2*np.ones(nruns, int);

OBCx = 1
OBCy = 1

PYROCHLORE = 0
if PYROCHLORE:
    OBCx = 1
    OBCy = 1

#EXCHANGE#
<<<<<<< HEAD
tl     = np.linspace(0.85, 1.05, nruns)#np.ones(nruns)
=======
tl     = np.array([0, 0.25, 0.5, 0.75, 1, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 2, 3, 4, 5, 6, 7, 10, 20, 100]) #np.linspace(0, 5, nruns)#np.ones(nruns)
>>>>>>> ef3c47272fadb627936e0457744360c245e4d47c
tr     = tl
Jzl    = -np.ones(nruns)
Jzr    = Jzl
Jpml   = np.zeros(nruns) #-(1*np.logspace(0, np.log10(2), nruns)-np.ones(nruns))
Jpmr   = Jpml

OBC = OBCx             #1 for open boundary conditions. 0 for periodic boundary conditions.
if OBCx and not OBCy:
    tr = 2*tr
    Jzr = 2*Jzr
    Jpmr = 2*Jpmr

if not OBCx and OBCy:
    tr = 0.5*tr
    Jzr = 0.5*Jzr
    Jpmr = 0.5*Jpmr

EIGVECS = 1         #Compute eigenvectors?
ZEROCORR = 1        #Compute all correlations for site zero?
NNCORR = 1          #Compute all nearest neighbour correlations?
MIDCORR = 1         #Compute all correlations for middle site?

RESETOLDFILES = 1



if mainproject:
    os.system("mkdir " + "Data/" + mainproject)
    totalproject = mainproject + "/" + project
else:
    totalproject = project

os.system("mkdir " + "Data/" + totalproject)

descfile = open("Data/" + totalproject + "/description.txt", 'w')
descfile.write(description)

runsub = open("Data/" + totalproject + "/run.sub", "w")

runsub.write("#!/bin/sh\n")
runsub.write("#SBATCH --account=nn4563k\n")
runsub.write("#SBATCH --time=" + time + "\n")
runsub.write("#SBATCH --partition=bigmem\n")
runsub.write("#SBATCH --mem-per-cpu=16999M\n")
runsub.write("#SBATCH --ntasks=1\n")
runsub.write("#SBATCH --signal=B:USR1@60\n")

if runsame:
    runsub.write("#SBATCH --array={}-{}%{}\n".format(runmin, runmax, runsame))
else:
    runsub.write("#SBATCH --array={}-{}\n".format(runmin, runmax))
runsub.write("#SBATCH --job-name="+jobname+"\n")
runsub.write("#SBATCH --error=stderr.dat\n")
runsub.write("#SBATCH --output=stdout.dat\n")
runsub.write("#set -o errexit\n")
runsub.write("\n")
runsub.write("\n")
runsub.write('#cleanup "rsync -av $SCRATCH/ $SUBMITDIR/ --exclude=stdout.dat --exclude=stderr.dat"\n')
runsub.write("#cd $SUBMITDIR\n")
runsub.write("#rsync -av $SUBMITDIR/ $SCRATCH/ --exclude=rundir\n")
runsub.write("#cd $SCRATCH\n")
runsub.write("echo Running program.....\n")
runsub.write("$HOME/Documents/OBCZRH/Code/program " + totalproject + " $SLURM_ARRAY_TASK_ID\n")





for run in range(runmin, nruns + runmin):
    #delta = np.ones(num)*(np.logspace(dmin, dmax, num)[run])
    if run < 10:
        os.system("mkdir " + "Data/" + totalproject + "/Run00" + str(run))
        outfile = open("Data/" + totalproject + "/Run00" + str(run) + "/parameters.txt", 'w')
    elif run < 100:
        os.system("mkdir " + "Data/" + totalproject + "/Run0" + str(run))
        outfile = open("Data/" + totalproject + "/Run0" + str(run) + "/parameters.txt", 'w')
    elif run < 1000:
        os.system("mkdir " + "Data/" + totalproject + "/Run" + str(run))
        outfile = open("Data/" + totalproject + "/Run" + str(run) + "/parameters.txt", 'w')
    else:
        print("generate project: Run number is bigger than 999")
        exit(1)

    outfile.write("Nsites = ")
    outfile.write(str(Nsites[run-runmin]))
    outfile.write("\n")
    outfile.write("Nh = ")
    outfile.write(str(Nh[run-runmin]))
    outfile.write("\n")
    outfile.write("\n")

    outfile.write("tl = ")
    outfile.write(str(tl[run-runmin]))
    outfile.write("\n")
    outfile.write("tr = ")
    outfile.write(str(tr[run-runmin]))
    outfile.write("\n")
    outfile.write("Jzl = ")
    outfile.write(str(Jzl[run-runmin]))
    outfile.write("\n")
    outfile.write("Jzr = ")
    outfile.write(str(Jzr[run-runmin]))
    outfile.write("\n")
    outfile.write("Jpml = ")
    outfile.write(str(Jpml[run-runmin]))
    outfile.write("\n")
    outfile.write("Jpmr = ")
    outfile.write(str(Jpmr[run-runmin]))
    outfile.write("\n")
    outfile.write("\n")

    outfile.write("OBC = ")
    outfile.write(str(OBC))
    outfile.write("\n")

    outfile.write("PYROCHLORE = ")
    outfile.write(str(PYROCHLORE))
    outfile.write("\n")

    outfile.write("EIGVECS = ")
    outfile.write(str(EIGVECS))
    outfile.write("\n")

    outfile.write("ZEROCORR = ")
    outfile.write(str(ZEROCORR))
    outfile.write("\n")

    outfile.write("NNCORR = ")
    outfile.write(str(NNCORR))
    outfile.write("\n")

    outfile.write("NNCORR = ")
    outfile.write(str(MIDCORR))
    outfile.write("\n")

    outfile.write("RESETOLDFILES = ")
    outfile.write(str(RESETOLDFILES))
    outfile.write("\n")

    #outfile.write(str(BOOST))


for i in range(runmin, runmax+1):
    if mainproject:
        shellfile = open('Jobs/' + mainproject + project + str(i) + '.sh', 'w')
    else:
        shellfile = open('Jobs/' + project + str(i) + '.sh', 'w')
    shellfile.write('#!/bin/bash\n')
    shellfile.write('nice -' + str(NICE) + ' ~/Documents/OBCZRH/Code/program ' + totalproject + ' ' + str(i) + ' &\n')
    if mainproject:
        shellfile.write('rm ' + '~/Documents/OBCZRH/Jobs/' + mainproject + project + str(i) + '.sh')
    else:
        shellfile.write('rm ' + '~/Documents/OBCZRH/Jobs/' + project + str(i) + '.sh')
    shellfile.close()
    if mainproject:
        os.system('chmod +x Jobs/' + mainproject + project + str(i) + '.sh')
    else:
        os.system('chmod +x Jobs/' + project + str(i) + '.sh')
