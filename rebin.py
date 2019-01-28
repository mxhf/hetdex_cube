# coding: utf-8
from __future__ import print_function

import pylab
from astropy.io import fits
import sys
from scipy import *
from scipy import optimize
import spectrum
import numpy
import optparse
import PolygonIntersect
from matplotlib import pyplot as plt
import spectrum
from astropy.io import ascii
from astropy.table import Table
import os
import glob
import numpy as np

from scipy.interpolate import UnivariateSpline

import srebin
import sys
from astropy.table import Column, vstack
from astropy.io import ascii
import argparse
import pickle
import glob

def get_rebinned(hdu, extensions, start = 3494.74, step =  1.9858398, stop = 5500.):
    #start,stop = 3503.9716796, 5396.477
    N = int( np.ceil( (stop - start)/step ) )

    rebinned = {}
    wl = hdu['wavelength'].data
    for ext in extensions:
        new = np.zeros([wl.shape[0], N])

        for i in range(wl.shape[0]):
            w = wl[i,:]
            f = hdu[ext].data[i,:]
            start = start
            step =  step
            stop = stop
            lw, lf = srebin.linlin(w, f, start, step, stop, dowarn=False)
 
            # hack as they may not all necessareyly have the same length
            new[i,:min(N, len(lf))] = lf[:min(N, len(lf))]

        rebinned[ext] = new
    return lw, rebinned


def read_shotlist(shotlist_file):
    return ascii.read(shotlist_file, format="fixed_width")


parser = argparse.ArgumentParser(description='Linearely rebin multifits files.')
parser.add_argument('--basepath', default="../reductions")
parser.add_argument('--force_rebin', action="store_true",
                            help='Force rebinning rather then using prior cached rebinning results.')
parser.add_argument('--shot', type=str, default="", 
                            help='Actual shots to rebin.')



args = parser.parse_args()

basepath = args.basepath
prefix = ""
extensions = ["spectrum", "sky_subtracted", "sky_spectrum", "fiber_to_fiber"]



night,shotid = args.shot.split("v")

pattern="{}/{}/virus/virus0000{}/exp??".format(basepath,night,shotid)
dd = glob.glob(pattern)
exposures = [ os.path.split(d)[1] for d in dd ]

print("Found {} exposures ...".format( exposures ) )
###############################################################################
# read spectra
###############################################################################

for exp in  exposures:
    shot =  "{}v{}".format(night,shotid)

    pattern="{}/{}/virus/virus0000{}/{}/virus/multi_???_???_???_??.fits".format(basepath,night,shotid, exp)

    ff = glob.glob(pattern)
    print(pattern)
    for filename  in ff:
        rebin_path = "{}v{}/{}".format(night, shotid, exp)
        __,t  = os.path.split( filename )
        rebin_filename = t.replace(".fits","_rebin.pickle")
        rebin_file_path = os.path.join(rebin_path,rebin_filename)
        if os.path.exists( rebin_file_path ) and not args.force_rebin:
            print("{} already exists, skipping ...".format(rebin_file_path))
            continue
        else:
            print("Read & rebin:", filename)
            hdu = fits.open(filename)
            lw, rebinned = get_rebinned(hdu, extensions, start = 3494.74, step =  1.9858398, stop = 5500.)
            try:
                os.makedirs(rebin_path)
            except:
                pass
            with open(os.path.join(rebin_path,rebin_filename), 'wb') as f:
                # Pickle the 'data' dictionary using the highest protocol available.
                pickle.dump((lw, rebinned), f, pickle.HIGHEST_PROTOCOL)

cmd = 'tar cvzf {}.tar.gz {}'.format(args.shot,args.shot)
os.system(cmd)
cmd = 'rm -rf {}'.format(args.shot)
os.system(cmd)

