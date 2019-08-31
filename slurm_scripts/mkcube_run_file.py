
import sys

cmd="python ../../src/hetdex_cube/cube.py --no_pca_sky --shiftsdir karls_shifts  --ifuslot {ifu} --shotlist ../shifts/repeat_{field}.txt -o outcube_{field}_{ifu}.fits.gz"

listfile = sys.argv[1]
field    = sys.argv[2]
runfile  = sys.argv[3]
with open(listfile,'r') as f:
    ll = f.readlines()

n = len(ll)
with open(runfile, 'w') as f:
    for l in ll:
        ifu = l.replace("\n","")
        f.write(cmd.format(ifu = ifu, field = field) + "\n")

    for i in range(n,100):
        f.write("echo 'Filler line.'\n")
