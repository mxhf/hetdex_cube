from astropy.io import fits
import spectrum
from matplotlib import pyplot as plt
import numpy as np
from astropy.stats import biweight_midvariance as bwtmv
import argparse
import sys

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--incube", type=str, 
                    help="Input data cube.")
parser.add_argument("-t", "--threshold", type=int, default=-0.005,
                    help="Upper threshold for data to be considered in the noise statistics.")
parser.add_argument("-o", "--outfile", type=str,
                    help="Input data cube.")


args = parser.parse_args(sys.argv[1:])

th = args.threshold
s = spectrum.readSpectrum(args.incube)
#ns = spectrum.readSpectrum("data/sf2ncoutcube_COSMOSB_022.fits.gz")



# compute noise in negative signal as unction of wavelength
dd = []
for sliceid in range(s.data.shape[0]):
    vv = s.data[sliceid][s.data[sliceid]!=0.].flatten()
    vv = np.append( vv[vv<th], -vv[vv<th] )
    bwstd = np.sqrt( bwtmv(vv) ) 
    std   =  np.std(vv) 
    dd.append([sliceid, std, bwstd])
dd = np.array(dd)

# fit model
p = np.polyfit(s.grid(), dd[:,2], deg=3)


# diagnostic plot
f = plt.figure(figsize=[10,5])

ax = plt.subplot(121)
sliceid = 500
vv = s.data[sliceid][s.data[sliceid]!=0.].flatten()

vv = np.append( vv[vv<th], -vv[vv<th] )
__ = plt.hist( vv , bins = np.arange(-.5,.5,.02), alpha=.5)
#__ = plt.hist( ns.data[sliceid].flatten() , bins = np.arange(-.5,.5,.02), alpha=.5)
#plt.text(.1,.9,"IFU {}".format(IFU), transform=ax.transAxes)

plt.subplot(122)
#plt.plot(dd[:,0],dd[:,1]*5.,'-')
plt.plot(s.grid(),dd[:,2],'-')
plt.plot(s.grid(),dd[:,2],'-')
plt.plot( s.grid(), np.polyval(p, s.grid()) )
plt.xlabel("slice")
plt.ylabel("std, bwstd")
plt.ylim([0.,.4])

f.tight_layout()

fnplot = args.incube.replace(".fits.gz","detect_noise_model.pdf")
plt.savefig( fnplot )
print("Wrote {}.".format(fnplot))

np.savetxt(args.outfile, p)

print("Wrote {}.".format(args.outfile))