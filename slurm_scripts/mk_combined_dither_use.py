
import sys
from astropy.table import Column, vstack
from astropy.io import ascii

shotlist_file = sys.argv[1]
out_file = sys.argv[2]

tables = []
with open(shotlist_file, 'r') as fin:
    ll = fin.readlines()
    for l in ll:
        if l.strip().startswith("#"):
            continue
        tt = l.split()
        night   = tt[0]
        shotid = tt[1]
        dithall_filename = "{}v{}/dithall.use".format(night, shotid)
        print("Reading {}".format(dithall_filename))
        t = ascii.read(dithall_filename)
        cn = Column(name="night", data=[night]*len(t) )
        cs = Column(name="shotid", data=[shotid]*len(t) )
        t.add_column(cs,index=0)
        t.add_column(cn,index=0)
        tables.append(t)


T = vstack(tables)
T.write(out_file, overwrite=True, format="ascii.fixed_width")
