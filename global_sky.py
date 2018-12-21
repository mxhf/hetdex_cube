
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
from astropy.table import Table, Column
from astropy import table
import os
import glob
import numpy as np
from scipy.interpolate import UnivariateSpline
import srebin
import sys
from astropy.table import Column, vstack
from astropy.io import ascii

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
    return Table(ascii.read(args.shotlist, format="fast_no_header"), names=["shots"])


def combine_dithall(shotlist, shifts_dir):
    tables = []
    for shot in shotlist["shots"]:
            dithall_filename = "{}/dithall.use".format(shot)
            filename = os.path.join(shifts_dir, dithall_filename)
            print("Reading {}".format( filename ))
            t = ascii.read(filename)
            night, shotid = shot.split("v")
            cn = Column(name="night", data=[night]*len(t) )
            cs = Column(name="shotid", data=[shotid]*len(t) )
            t.add_column(cs,index=0)
            t.add_column(cn,index=0)
            tables.append(t)

    T = vstack(tables)
    return T




RA0 = None 
DEC0 = None
FIBERD = 1.5
nx = None
ny = None
pixelsize = .5

fiberA = pi*(FIBERD/2.)**2.

import argparse

parser = argparse.ArgumentParser(description='Build a hetdex cube.')
#parser.add_argument('--basepath', default="/work/03946/hetdex/maverick/red1/reductions")
parser.add_argument('--basepath', default="../reductions")
parser.add_argument('--shiftsdir', type=str, default="../shifts",
                            help='Directory that contains the astrometric solutions for each shot.')
parser.add_argument('--dither_use', type=str,
                            help='Combined dithall use file.')
parser.add_argument('--force_rebin', action="store_true",
                            help='Force rebinning rather then using prior cached rebinning results.')
parser.add_argument('--shotlist', type=str,
                            help='List of actual shots to use.')
args = parser.parse_args()

basepath = args.basepath
prefix = ""
extensions = ["sky_subtracted", "sky_spectrum", "fiber_to_fiber"]

shotlist = read_shotlist(args.shotlist)
t = combine_dithall(shotlist, args.shiftsdir)

# read dithall.use
#filebase = {}
spectra = {}
#exposure_times = []
#allspec = []
per_amp_sky_spectra = [] # holds for each fiber spectrum
                 # the corresponing amplifier wide sky (median accorss all fibers after correcting for fiber_to_fiber)
wlgrid = None
names = ["count", "amplifier", "fiberid", "ra", "dec", "shot", "night", "shotid", "exp"]
dtype = [int, 'U2', int, float, float, 'U12', 'U8', 'U3', 'U6']
fibers = Table(names=names, dtype=dtype)

count = 0
fid = -1

import pickle

###############################################################################
# read spectra
###############################################################################


def amp(x):
    return x[18:20]

camp = map(amp, t["multifits"]) 
t.add_column(Column(camp, name='amp') )


ut1 = table.unique(t, keys=['night', 'shotid', 'exposure', 'ifuslot', 'amp'])
for r in ut1:
    mf = prefix + r["multifits"]
    #_ifuslot = r["ifuslot"].replace("ifu","")
    #tt = mf.split("_")
    #fiberid = int( tt[5][:3] ) - 1
    #amplifier = mf[-10:-8]
    exp = r["exposure"]
    night = r["night"]
    shotid = "{}".format( r["shotid"] )
    shot =  "{}v{}".format(night,shotid)
    if (not (shot in shotlist["shots"])):
        continue
    filename = mf[:-8] + ".fits"
    #id = int(tt[1])
    #x = r["ra"]
    #y = r["dec"]
    #date = r["timestamp"][:8]
    filename = mf[:-8] + ".fits"
    s= "{}/virus/virus0000{}/{}/virus/{}".format(night, shotid, exp, filename)
    #ff.append(filename)
    path = os.path.join( basepath, s )

    # read spectrum if it was not read before
    if not path in spectra:
        rebin_path = "rebin/{}v{}/{}".format(night, shotid, exp)
        rebin_filename = filename.replace(".fits","_rebin.pickle")
        rebin_file_path = os.path.join(rebin_path,rebin_filename)
        if os.path.exists( rebin_file_path ) and not args.force_rebin:
            # already rebinned?
            with open( rebin_file_path , 'rb') as f:
                print("Found previously rebinned data {}".format( rebin_file_path )) 
                # The protocol version used is detected automatically, so we do not
                # have to specify it.
                lw, rebinned = pickle.load(f)
                wlgrid = lw
        else:
            print("Read & rebin:", path)
            if not os.path.exists( path ):
                print("WARNING: Did not find {}".format(path))
            else:
                hdu = fits.open(path)
                lw, rebinned = get_rebinned(hdu, extensions, start = 3494.74, step =  1.9858398, stop = 5500.)
                if type(wlgrid) == type(None):
                    wlgrid = lw
                try:
                    os.makedirs(rebin_path)
                except:
                    pass
                with open(os.path.join(rebin_path,rebin_filename), 'wb') as f:
                    # Pickle the 'data' dictionary using the highest protocol available.
                    pickle.dump((lw, rebinned), f, pickle.HIGHEST_PROTOCOL)

    per_amp_sky_spectra.append( np.nanmedian( rebinned['sky_spectrum']/rebinned['fiber_to_fiber'], axis=0 ) )

per_amp_sky_spectra = np.array(per_amp_sky_spectra)
per_shot_sky_spectra = {}
# now average all sky spectra for each shot
ut2 = table.unique(t, keys=['night', 'shotid'])
for r in ut2:
    ii = ( ut1['night'] == r['night']) * (ut1['shotid'] == r['shotid'])
    ff = np.nanmedian(sky_spectra[ii],axis=0)
    sout = Table([lw, ff], names=['wavelength', 'counts'], dtype=[float,float])
    sout.write("{}v{}_sky.fits".format(r['night'], r['shotid']), format="fits")
    
