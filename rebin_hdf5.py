#!/usr/bin/env python

from __future__ import print_function

import sys
import os

import numpy as np
import argparse
import pickle
#from astropy.table import Table
from scipy.interpolate import UnivariateSpline
import tables as tb

import srebin


def get_rebinned(ww, ff, start = 3494.74, step =  1.9858398, stop = 5500.):
    N = int( np.ceil( (stop - start)/step ) )
    _lw, _lf = srebin.linlin(ww, ff, start, step, stop, dowarn=False)
    
    lw = (np.arange(N) * step) + start
    lf = np.zeros_like(lw)
    
    lf[:min(N, len(_lf))] = _lf[:min(N, len(_lf))]

    if not _lw[0] == lw[0]:
        print("Binning is broken, please check.")
        return np.nan
    
    return lw, lf

def rebin_IFU(t,  expnum, ifuslot, rebin_path, force_rebin=False):
    print("Processing expnum {} ifuslot {} ".format(expnum, ifuslot.decode("utf-8")) )
    amplifiers = \
        np.unique(np.unique( [ x['amp'] for x in t.where("""(expnum == {expnum}) & (ifuslot == {ifuslot})""".format(expnum=expnum, ifuslot=ifuslot))] ))

    for amp in amplifiers:

        rebin_filename = "multi_xxx_{ifuslot}_xxx_{amp}_rebin.pickle".format(ifuslot=ifuslot.decode("utf-8"),amp=amp.decode("utf-8"))
        rebin_file_path = os.path.join(rebin_path,rebin_filename)

        if os.path.exists( rebin_file_path ) and not force_rebin:
            print("{} already exists, skipping ...".format(rebin_file_path))
            continue
        else:    
            rebinned = {}

            idx = [ x['fibidx'] for x in t.where( """(expnum == {expnum}) & (ifuslot == {ifuslot}) & (amp == {amp})""".format(expnum=expnum, ifuslot=ifuslot, amp=amp) ) ]
            www = [ x["wavelength"] for x in t.where( """(expnum == {expnum}) & (ifuslot == {ifuslot}) & (amp == {amp})""".format(expnum=expnum, ifuslot=ifuslot, amp=amp) ) ]

            for ext in extensions:
                fff = [ x[ext] for x in t.where( """(expnum == {expnum}) & (ifuslot == {ifuslot}) & (amp == {amp})""".format(expnum=expnum, ifuslot=ifuslot, amp=amp) ) ]
                data = []
                print("Rebinning expnum {} , ifuslot {} , amp {} , ext {} ".format( expnum, ifuslot.decode("utf-8"), amp.decode("utf-8"), ext) )

                for i,ww,ff in zip(idx,www,fff):

                    lw, lf = get_rebinned(ww, ff, start = 3494.74, step =  1.9858398, stop = 5500.)
                    data.append(lf)

                rebinned[ext] = np.array(data)

            try:
                os.makedirs(rebin_path)
            except:
                pass

            with open(rebin_file_path, 'wb') as f:
                # Pickle the 'data' dictionary using the highest protocol available.
                print("Writing ", rebin_file_path)
                pickle.dump((lw, rebinned), f, pickle.HIGHEST_PROTOCOL)
                #sys.exit(1)
    return True


# parse arguments
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-s', '--shot', default='20180124v010', required=True)
parser.add_argument('-i', '--ifuslot', default='')
parser.add_argument('--data_dir', default="/work/03946/hetdex/hdr1/reduction/data/")
parser.add_argument('--force_rebin', action='store_true')
parser.add_argument('--output_dir', default='/work/04287/mxhf/maverick/hetdex/rebin2')

args = parser.parse_args()


# extensions that will be rebinned
extensions = ["spectrum", "sky_subtracted", "sky_spectrum", "fiber_to_fiber"]
filename = os.path.join(args.data_dir, '{}.h5'.format(args.shot))

# load shot
fileh = tb.open_file(filename)
t = fileh.root.Data.Fibers

# exposures  to operate on
exposures = [1,2,3]       

# find which ifuslots we have
ifuslots = np.unique(np.unique(  [x['ifuslot'] for x in t]  ))

if args.ifuslot != '':
	ifuslots = [ args.ifuslot.encode() ]

#execute
for ifuslot in ifuslots:
        for expnum in exposures:
                rebin_path = os.path.join(args.output_dir, "{}/exp{:02d}".format(args.shot, int(expnum)) )
                rebin_IFU(t, expnum, ifuslot, rebin_path, args.force_rebin)
