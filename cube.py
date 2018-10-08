
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
import srebin

# In[146]:


def circle(x,y,size):
        tt = arange(0.,2.*pi, pi/18.)
        xx = size/2. * cos(tt) + x
        yy = size/2. * sin(tt) + y
        return xx,yy

def pixel(x,y,size):
    xx,yy=[],[]
    xx.append(x-size/2.)
    yy.append(y-size/2.)

    xx.append(x+size/2.)
    yy.append(y-size/2.)

    xx.append(x+size/2.)
    yy.append(y+size/2.)

    xx.append(x-size/2.)
    yy.append(y+size/2.)

    return xx,yy

def create_3D_header(xc, yc, rad, decd, pixelsize,start,step):
    h = fits.Header()
    #h.update("WCSDIM  ", 3, "");
    #h.update("WAT0_001", "system=image", "");
    h["CTYPE3"] = "Wave"
    h["CRPIX3"] = 1.
    h["CRVAL3"] = start
    h["CDELT3"] = step
    #h["WAT1_001"] = wtype=tan
    h["CTYPE1"] = "RA--TAN"
    h["CRPIX1"] = xc
    h["CRVAL1"] = rad
    h["CDELT1"] = - pixelsize/3600.
    h["CUNIT1"] = "deg"
    #h["WAT2_001"] = wtype=tan
    h["CTYPE2"] = "DEC--TAN"
    h["CRPIX2"] = yc
    h["CRVAL2"] = decd
    h["CDELT2"] = pixelsize/3600.
    h["CUNIT2"] = "deg"
    return h

def create_2D_header(xc, yc, rad, decd, pixelsize):
    h = fits.Header()
    #h.update("WAT1_001", "wtype=tan axtype=ra", "");
    h["CTYPE1"] = "RA--TAN"
    h["CRPIX1"] = xc
    h["CRVAL1"] = "rad"
    h["CDELT1"] = - pixelsize/3600.
    h["CUNIT1"] = "deg"
    #h["WAT2_001"] = wtype=tan
    h["CTYPE2"] = "DEC--TAN"
    h["CRPIX2"] = yc
    h["CRVAL2"] = decd
    h["CDELT2"] = pixelsize/3600.
    h["CUNIT2"] = "deg"
    return h

def create_I_header():
    h = fits.Header()
    #h.update("WAT1_001", "wtype=tan axtype=ra", "");
    h["CTYPE1"] = "pixel"
    h["CRPIX1"] = 1
    h["CRVAL1"] = 1.
    h["CDELT1"] = 1.
    h["CUNIT1"] = "px"
    #h["WAT2_001"] = wtype=tan
    h["CTYPE2"] = "fiber"
    h["CRPIX2"] = 1
    h["CRVAL2"] = 1.
    h["CDELT2"] = 1.
    h["CUNIT2"] = "fib"
    return h

def tan_dir_sci(RA0, DEC0, PA0, RA,DEC, quiet=True):
    """
    Calculates pixel positions in the IFU for a given set of RA and DEC coordinates.

    Input

    IFU_RA0, IFU_DEC0 = IFU zero coordinates (in deg)
    RA, DEC = 1D arrays of RA and DEC coordinates
    astrom_terms = parameters for the astrometric solution
    quiet = Boolean, disable output (standard: False)

    Returns:
    pixx, pixy = two 1D arrays containing x and y coordinates (in deg!) for the given RAs and DECs
    """
    if not quiet: print("IFU_RA0, IFU_DEC0 = {},{}".format( IFU_RA0, IFU_DEC0) )

    rRA0  = RA0*pi/180.
    rDEC0  = DEC0*pi/180.

    rPA0 = PA0*pi/180.

    if not quiet: print("[tan_dir_sci] IFU_RA0, IFU_DEC0 = {},{}".format( IFU_RA0, IFU_DEC0 ) )

    rRA = RA*pi/180.

    rDEC = DEC*pi/180.

    if not quiet: print("[tan_dir_gui] rRA, rDEC: {}, {}".format( rRA, rDEC ) )

    # eq 9 in Greisen AIPSMEMO27
    L = cos(rDEC)*sin(rRA - rRA0)/        (sin(rDEC)*sin(rDEC0) + cos(rDEC)*cos(rDEC0)*cos(rRA -         rRA0))

    M = (sin(rDEC)*cos(rDEC0) -         cos(rDEC)*sin(rDEC0)*cos(rRA - rRA0))/(sin(rDEC)*        sin(rDEC0) + cos(rDEC)*cos(rDEC0)*cos(rRA - rRA0))

    # eq 5 in Greisen AIPSMEMO27
    pixx = L*cos(rPA0) + M*sin(rPA0)
    pixy = M*cos(rPA0) - L*sin(rPA0)

    return - pixx/pi*180.*3600., pixy/pi*180.*3600.

def findZeroPixRaDec(x,y, RA0, DEC0):
    """
    Find RA and DEC coordinates that correspond to the given
    pixel coordinates.
    We do this here the cheap way. We use a nonlinear fit 
    rather than working out the inverse transformation.
    """
    def peval(p,RA0, DEC0):
        RA,DEC = p
        return tan_dir_sci(RA0, DEC0, 0., RA, DEC, quiet=True)

    def resid(p,x,y,RA0, DEC0):
        xt,yt = peval(p,RA0, DEC0)
        return (xt-x), (yt-y)

    p0 = [RA0,DEC0]
    bestfit = optimize.leastsq(resid, p0, args=(x,y,RA0,DEC0) )
    return bestfit[0]

def hms2deg(hms):
    tt = hms.split(":")
    h,m,s = float(tt[0]), float(tt[1]), float(tt[2])
    return h*15. + m/4. + s/240.

def dms2deg(dms):
    tt = dms.split(":")
    d,m,s = float(tt[0]), float(tt[1]), float(tt[2])
    return d + m/60. + s/3600.


# In[147]:


def get_rebinned(hdu, extensions=['spectrum'], start = 3494.74, step =  1.9858398, stop = 5500.):
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


# In[148]:



# In[149]:

FIBERD = 1.5
fiberA = pi*(FIBERD/2.)**2.

RA0 = None 
DEC0 = None
pa = 264.116951
FIBERD = 1.
nx = None
ny = None
pixelsize = 1.


import argparse

parser = argparse.ArgumentParser(description='Build a hetdex cube.')
parser.add_argument('--basepath', default="/work/03946/hetdex/maverick/red1/reductions")
parser.add_argument('--pa', type=float, default=0.,
                            help='Position angle for cube.')
parser.add_argument('--shots_list', type=str, default="",
                            help='List of actual shots to use.')
parser.add_argument('dither_use', type=str,
                            help='Dither.use file')
parser.add_argument('ifuslot', type=str, default = "022",
        help='IFUslot to create cube for. ')

args = parser.parse_args()

pa=args.pa
fiberpos = args.dither_use
#"/work/04287/mxhf/maverick/sci/panacea/shifts/deep_33.182404_262.562869.use"
ifuslot = args.ifuslot
#"022"
basepath = args.basepath
#"/work/03946/hetdex/maverick/red1/reductions"
prefix = ""
extension = "sky_subtracted"


shots = []
if args.shots_list != "":
    with open(args.shots_list, 'r') as f:
        ll = f.readlines()
    for l in ll:
        if l.strip().startswith("#") or l.strip() == "":
            continue
        shots.append(l.strip())
# In[150]:





# In[151]:


# read dithall.use
filebase = {}
spectra = {}
exposure_times = []
allspec = []

wlgrid = None
fibers = Table(names=["count","ra", "dec"])

count = 0
fid = -1

print("Reading {}".format(fiberpos))
t = ascii.read(fiberpos,format="fixed_width")

import pickle

for r in t:
    mf = prefix + r["multifits"]
    _ifuslot = r["ifuslot"].replace("ifu","")
    if not ifuslot == _ifuslot:
        continue
    tt = mf.split("_")
    fiberid = int( tt[5][:3] ) - 1
    exp = r["exposure"]
    night = r["night"]
    shotid = "{:03d}".format( r["shotid"] )
    shot =  "{}v{}".format(night,shotid)
    if len(shots) > 1  and not shot in shots:
        print("Skipping {} as it is not listed in shot list (--shots_list).".format(shot))
        continue
    filename = mf[:-8] + ".fits"
    id = int(tt[1])
    x = r["ra"]
    y = r["dec"]
    date = r["timestamp"][:8]
    filename = mf[:-8] + ".fits"
    s= "{}/virus/virus0000{}/{}/virus/{}".format(night, shotid, exp, filename)
    path = os.path.join( basepath, s )
    # read spectrum if it was not read before
    if not path in spectra:
        rebin_path = "rebin/{}v{}/{}".format(night, shotid, exp)
        rebin_filename = filename.replace(".fits","_rebin.pickle")
        rebin_file_path = os.path.join(rebin_path,rebin_filename)
        if os.path.exists( rebin_file_path ):
            # already rebinned?
            with open( rebin_file_path , 'rb') as f:
                print("Found previously rebinned data {}".format( rebin_file_path )) 
                # The protocol version used is detected automatically, so we do not
                # have to specify it.
                lw, rebinned = pickle.load(f)
                wlgrid = lw
        else:
            print("Read & rebin:", path)
            hdu = fits.open(path)
            lw, rebinned = get_rebinned(hdu, extensions=[extension], start = 3494.74, step =  1.9858398, stop = 5500.)
            if type(wlgrid) == type(None):
                wlgrid = lw
            try:
                os.makedirs(rebin_path)
            except:
                pass
            with open(os.path.join(rebin_path,rebin_filename), 'wb') as f:
                # Pickle the 'data' dictionary using the highest protocol available.
                pickle.dump((lw, rebinned), f, pickle.HIGHEST_PROTOCOL)

        spectra[path] = rebinned

    fibers.add_row([count,x,y])
    allspec.append(spectra[path][extension][fiberid,:])
    count += 1





# In[152]:


###############################################################################
# determine a pixel grid
###############################################################################

if RA0 == None or DEC0 == None:
    RA0,DEC0 = numpy.mean(fibers["ra"]), numpy.mean(fibers["dec"])
print("Set tangent point for projection to RA0 = %.6f and DEC0  = %.6f." % (RA0, DEC0)) 

fxx,fyy = tan_dir_sci(RA0, DEC0, pa, fibers["ra"], fibers["dec"], quiet=True)


maxx = max(fxx)+FIBERD/2.
maxy = max(fyy)+FIBERD/2.
minx = min(fxx)-FIBERD/2.
miny = min(fyy)-FIBERD/2.

if nx == None or ny == None:
    # Automatically determine pixel gridsize
    nx = int( round( ( maxx - minx ) / pixelsize ))
    ny = int( round( ( maxy - miny ) / pixelsize ))
    # create list of all pixel center coordinates 
    xx=arange(nx)*pixelsize + minx + pixelsize/2.
    yy=arange(ny)*pixelsize + miny + pixelsize/2.
else:
    # use user-defined grid insted
    xx=arange(nx)-(nx-1.)/2.
    yy=arange(ny)-(ny-1.)/2.
    xx *= pixelsize
    yy *= pixelsize

X,Y=meshgrid(xx,yy)
pixels = zeros( [ len(X.flatten()) ,3]  ) 
pixels[:,0] = arange(len(pixels))
pixels[:,1] = X.flatten()
pixels[:,2] = Y.flatten()

PLOT = False
if PLOT:
    # plotting
    s = plt.subplot() 
    for p in pixels:
        xx,yy = pixel(p[1],p[2], pixelsize)
        s.fill( xx, yy, facecolor='none',linewidth=0.2 )

    for f in zip( fxx,fyy ):
        xx,yy = circle(f[0],f[1], FIBERD)
        s.fill( xx, yy, facecolor='none',linewidth=0.2,edgecolor='blue' )

    s.set_xlabel("x (\")")
    s.set_ylabel("y (\")")
    s.axis('equal')
    


# In[153]:


###############################################################################
# calculate pixel/fiber intersections
###############################################################################
print( "Calculating fiber/pixel weight maps..." )
I = numpy.zeros([len(fibers),len(pixels)], dtype=float) # I is the pixel fiber intersection matrix

fib_range = numpy.arange( len(fibers) )
pix_range = range( len(pixels) )

NNN = 0

for ip in pix_range:
    p = pixels[ip]

    px,py = p[1],p[2]
    #calculate distances of all fibers to the current pixel
    dd_sq = (fxx-px)**2. + (fyy-py)**2.    

    #find which fibers could possibly intersect with the current pixel
    # in the intersection filter, only fibers which overlap the pixel are considered.
    # We only look at fibers wich are not further than 
    # sqrt(2) * pixelsize/2 + fiberd/2 
    ii = ( dd_sq < (pixelsize/2. *  1.414 + FIBERD/2.)**2.)
    # create a polygon describing the current pixel
    ppxx,ppyy = pixel(px,py,pixelsize)
    pixel_poly = list( zip(ppxx,ppyy) )

    if any(ii):
        for ifib in fib_range[ii]:
            fx = fxx[ifib]
            fy = fyy[ifib]
            if True:
                NNN += 1
                fpxx,fpyy = circle(fx,fy,FIBERD)
                fiber_poly = list( zip(fpxx,fpyy) )
                fiber_array = PolygonIntersect.toPointsArray(fiber_poly)
                pixel_array = PolygonIntersect.toPointsArray(pixel_poly)
                # calculate intersection area
                iA = PolygonIntersect.intersectionArea(fiber_array, pixel_array)
                # Now, the flux of a given fiber (at a given wavelength)
                # will be assigned to a pixel weighted by the fraction of the 
                # total fiber area that is overlapping with the pixel.
                I[ifib,ip]  = iA#/fiberA


plot_count = 0.
for ip in pix_range:
    #PLOT = (ip > 395) and (ip < 400)
    PLOT = False
    if PLOT:
        print("pixel ", ip)
        s = pylab.axes() 
        p = pixels[ip]
        px,py = p[1],p[2]
        s.plot([px],[py],'s')
        ppxx,ppyy = pixel(px,py,pixelsize)
        pixel_poly = zip(ppxx,ppyy)
        s.fill( ppxx,ppyy, facecolor='none',edgecolor='k',linewidth=0.2 )

        #if any(ii):
        #   for j in jj[ii]:
        for ifib in fib_range:
            if I[ifib,ip] > 0.:
                fx = fxx[ifib]
                fy = fyy[ifib]

                fpxx,fpyy = circle(fx,fy,FIBERD)

                s.plot([fx],[fy],'.')
                s.fill( fpxx,fpyy, facecolor='none',edgecolor='b',linewidth=0.2 )
                s.text(fx,fy,'%.2f' % I[ifib,ip])

        #s.text(px,py,'%.2f' % W[i])
        pylab.show()

print("Done.")


# In[154]:


nI = I/fiberA
IT = nI.transpose()

wlstart, wlstop = wlgrid[0], wlgrid[-1]
allspec = np.array(allspec)



# In[155]:


###############################################################################
# Create Cube
###############################################################################
print("Creating cube...")
# with the help of the intersection matix
# the pixel values of each wavelength slice are simply the dot product of that matrix (transposed) in the vector
# of all the fiber values at a give wavelength.
cube = 1./np.zeros( [allspec.shape[1],ny,nx] )

W = sum(I,axis=0)/(pixelsize**2.)

# calculate median flux level of all fibers
kk = (wlgrid > wlstart) * (wlgrid < wlstop) 
mm = np.median(allspec[:,kk],axis=1)

nallspec = (allspec.transpose()/mm).transpose()


def stats(a):
    print("min,max = ", a.min, a.max() )
    print("mean = ", mean(a) )
    print("std = ", std(a) )

if True:
    #for iwl in range(allspec.shape[1])[361:366]: # for each wavelength
    for iwl in range(allspec.shape[1]): # for each wavelength

        if iwl % 100 == 0:
            s = "WL: %6d %14.6f" % ( iwl,wlgrid[iwl] )
            sys.stdout.write( '\r'* len(s))
            sys.stdout.write(s)
            sys.stdout.flush()

        if True:
            imf2 = IT.dot(allspec[:,iwl])
            ii = W > 0. # prevent zero div error
            imf2[ii] = imf2[ii]/W[ii]
            im2 = imf2.reshape(X.shape)
        else:
            #pp = zeros( len(pixels) )
            Icr = nI.copy()
            jj = allcrmasks[:,iwl] > 0. # masked fibers (for current wavelength) 

            for ip in pix_range:
                ii = nI[:,ip] > 0.   # intersecting fiber apertures with this pixel
                Icr[ii*jj,ip] = 0.
                #for ifib in fib_range[ii*jj]:  
                #   print "1 wlind %5d, wl %8.5f, pixel %6d" % (iwl, wlgrid[iwl], ip), "Rejecting fiber", ifib
                #   Icr[ifib,ip] = 0.   # disable outlier fibers 

                # ###### PLOTTING ####
                #PLOT = (ip > 395) and (ip < 400)
                PLOT = False 
                if PLOT:
                    s = pylab.axes() 
                    vmin=0
                    vmax=500
                    cmap = pylab.cm.jet
                    for ifib in fib_range:
                        if nI[ifib,ip] > 0.:
                            fx = fxx[ifib]
                            fy = fyy[ifib]

                            fpxx,fpyy = circle(fx,fy,FIBERD)

                            s.plot([fx],[fy],'.')
                            val = allspec[ifib,iwl]
                            c = (val-vmin)/(vmax-vmin)
                            s.fill( fpxx,fpyy, facecolor=cmap(c),edgecolor='b',linewidth=1. )
                            s.text(fx,fy,'%.2f' % nI[ifib,ip])
                    
                            if allcrmasks[ifib,iwl] > 0.: # outlier fibers
                                print( "2 wlind %5d, wl %8.2f, pixel %6d" % (iwl, wlgrid[iwl], ip), "Rejecting fiber", ifib )
                                s.text(fx,fy,'X',fontsize=20)

                    p = pixels[ip]
                    px,py = p[1],p[2]
                    s.plot([px],[py],'s')
                    ppxx,ppyy = pixel(px,py,pixelsize)
                    pixel_poly = zip(ppxx,ppyy)
                    s.fill( ppxx,ppyy, facecolor='none',edgecolor='k',linewidth=0.2 )

                    #s.text(px,py,'%.2f' % W[i])
                    pylab.show()
                # ###### END PLOTTING ####

            # normalize by respective areas
            #W = sum(Icr,axis=0)/(pixelsize**2.)
            nIcr=Icr/fiberA
            ITcr = nIcr.transpose()

            # add to all pixel the fiber fluxes qeighted by ther fractional overlap with a respective pixel
            im2crrej = ITcr.dot(allspec[:,iwl])

            # correct for the pixel coverage
            ii = W > 0. # prevent zero div error
            im2crrej[ii] = im2crrej[ii]/W[ii]
        
            # correct for "missing" flux due to cr rejection    
            crnorm = (ITcr.dot(mm)/IT.dot(mm))
            ii = crnorm > 0.
            im2crrej[ii] = im2crrej[ii]/crnorm[ii]

            im2crrej = im2crrej.reshape(X.shape)
            im2 = im2crrej

        cube[iwl] = im2



# In[156]:


wstart = wlgrid[0]
wstep  = wlgrid[1]-wlgrid[0]


###############################################################################
# save output
###############################################################################

#xc,yc are the cube pixel indices that
#correspond to RA0 and DEC0 and x = 0" and y = 0"
xc = -xx[0]/pixelsize + 1 
yc = -yy[0]/pixelsize + 1 

h2d = create_2D_header(xc, yc, RA0, DEC0, pixelsize)
hdu = fits.PrimaryHDU(W.reshape(X.shape),h2d)
hdu.writeto("pixel_weights.fits.zg",overwrite=True)

h = create_I_header()
hdu = fits.PrimaryHDU(I,h)
hdu.writeto("fiber_weights.fits.gz",overwrite=True)

h = create_3D_header(xc, yc, RA0, DEC0, pixelsize,wstart,wstep)
#if options.normexptime:
#    h.add_history("Cube data has been normalized by its exposure time.")
hdu = fits.PrimaryHDU(cube,h)
hdu.writeto("outcube_{}.fits.gz".format(ifuslot),overwrite=True)
