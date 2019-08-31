
import os
import sys
import numpy as np


cmd="""python /home1/04287/mxhf/hetdex/xpca/pca_sky6.py --sky_kappa={sky_kappa} --ncomp={ncomp} -i {ifu} -a LL LU -b LU --shotlist_pca {shotlist_pca} --shotlist_skyrecon {shotlist_skyrecon} --dir_rebin {dir_rebin}
python /home1/04287/mxhf/hetdex/xpca/pca_sky6.py --sky_kappa={sky_kappa} --ncomp={ncomp} -i {ifu} -a LL LU -b LL --shotlist_pca {shotlist_pca} --shotlist_skyrecon {shotlist_skyrecon} --dir_rebin {dir_rebin}
python /home1/04287/mxhf/hetdex/xpca/pca_sky6.py --sky_kappa={sky_kappa} --ncomp={ncomp} -i {ifu} -a RL RU -b RU --shotlist_pca {shotlist_pca} --shotlist_skyrecon {shotlist_skyrecon} --dir_rebin {dir_rebin}
python /home1/04287/mxhf/hetdex/xpca/pca_sky6.py --sky_kappa={sky_kappa} --ncomp={ncomp} -i {ifu} -a RL RU -b RU --shotlist_pca {shotlist_pca} --shotlist_skyrecon {shotlist_skyrecon} --dir_rebin {dir_rebin}"""

slurmtempl="run_pca.slurm"
ifulistfile = sys.argv[1]
dir_rebin = sys.argv[2]
shotlist_pca = sys.argv[3]
shotlist_skyrecon = sys.argv[4]
sky_kappa = 1.5
ncomp = 25

ifus = []
def loadlist(listfile):
    li = []
    with open(listfile,'r') as f:
        ll = f.readlines()
    for l in ll:
        li.append(l.strip())
    return li

ifus = loadlist(ifulistfile)

worklist = []
for i in ifus:
    worklist.append(i)

N = len(worklist)
ii = np.arange(N)
iii = np.array_split(ii,int(48))

print("Distributing {} ifus over {} jobs.".format(N,len(iii)))
for j,ii in enumerate(iii):
    runfile="run_pca_{:04d}.run".format(j)
    with open(slurmtempl, 'r') as f:
        s = f.read()
    s = s.replace("@@@TASKS@@@", str(len(ii)*4) )
    s = s.replace("@@@CONTROL_FILE@@@", os.path.join('slurms_pca',runfile))
    with open(os.path.join('slurms_pca',slurmtempl.replace(".slurm", "_{:04d}.slurm".format(j))),'w') as f:
        f.write(s)
    path = os.path.join('slurms_pca',runfile)
    with open(path, 'w') as f:
        for i in ii:
            ifu = worklist[i]
            f.write(cmd.format(ifu = ifu, dir_rebin = dir_rebin, shotlist_pca = shotlist_pca, shotlist_skyrecon = shotlist_skyrecon, sky_kappa = sky_kappa, ncomp = ncomp) + "\n")
        print("Writing {}".format( path ) )

