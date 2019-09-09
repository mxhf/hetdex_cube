
import os
import sys
import numpy as np

cmd="python3 ../../src/hetdex_cube/check_rebin.py  -s {shot}"
slurmtempl="run_check_rebin_hdf5.slurm"
shotlistfile = "shotlist_PCA_COSMOSABCD.txt"

shots = []
def loadlist(listfile):
    li = []
    with open(listfile,'r') as f:
        ll = f.readlines()
    for l in ll:
        li.append(l.strip())
    return li

shots = loadlist(shotlistfile)

worklist = []
for s in shots:
    worklist.append(s)

N = len(worklist)
ii = np.arange(N)
iii = np.array_split(ii,int(50))

print("Distributing {} tasks over {} jobs.".format(N,len(iii)))
for j,ii in enumerate(iii):
    runfile="run_check_rebin_hdf5_{:04d}.run".format(j)
    with open(slurmtempl, 'r') as f:
        s = f.read()
    s = s.replace("@@@TASKS@@@", str(len(ii)) )
    s = s.replace("@@@CONTROL_FILE@@@", os.path.join('slurms_check_rebin',runfile))
    with open(os.path.join('slurms_check_rebin',slurmtempl.replace(".slurm", "_{:04d}.slurm".format(j))),'w') as f:
        f.write(s)
    with open(os.path.join('slurms_check_rebin',runfile), 'w') as f:
        for i in ii:
            shot = worklist[i]
            f.write(cmd.format(shot = shot) + "\n")
        print("Writing {}".format( runfile) )

