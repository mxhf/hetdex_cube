from astropy.io import ascii
import sys

fninfile  = sys.argv[1]
fnshotlist = sys.argv[2]
fnoutfile = sys.argv[3]

t       = ascii.read(fninfile, format="fixed_width")
shotlist = ascii.read(fnshotlist)["shotid"]

ii = ["{}v{:03d}".format(r['night'],r['shot']) in shotlist for r in t]

tout = t[ii]
tout.write(fnoutfile,overwrite=True, format="ascii.fixed_width")
