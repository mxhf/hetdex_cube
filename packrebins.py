
from numpy import loadtxt
import glob
import os

ifus = ["025"]
fnshotlist =  "shotlist_test.txt"
dir_rebin="/scratch/04287/mxhf/rebin2"
packdir = "./pack_test_025"

with open(fnshotlist) as f:
    ll = f.readlines()

shotlist = [l.strip() for l in ll]

for shot in shotlist:
    for ifu in ifus:
        for exp in ["exp01","exp02","exp03"]:
            tdir = "{}/{}/{}".format(packdir, shot, exp)
            cmd = "mkdir -p {}".format(tdir)
            os.system(cmd)
            pattern = "{}/{}/{}/multi_xxx_{}_xxx_??_rebin.pickle".format(dir_rebin, shot, exp, ifu)
            ff = glob.glob(pattern)
            for f in ff:
                cmd = "cp {} {}/.".format(f, tdir)
                os.system(cmd)


cmd = "tar -cvf {}.tar.gz {}".format(packdir, packdir)
os.system(cmd)
print("Produced {}.tar.gz".format(packdir))
