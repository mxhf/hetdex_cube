#!/bin/bash
# automaticallly resubmit jobs that failed

# find all .o files for jobs that got killed
for f in `ls HETDEX.o*`; do echo $f `tail -n1 $f`; done | grep "Killed by signal 15" | awk '{print $1}' > killed15

# find corresponding run files
for k in `cat killed15`; do cat $k | grep -o "rebin_hdf5.py -s ............" *.run | awk '{print $3}'; done | uniq > killed152

# find corresponding slurm files
for k in `cat killed152`; do grep $k *.run; done | sed 's/:/ /g' | awk '{print $1}' | uniq | sed 's/\.run/\.slurm/g' > killed15_3

# resubmit
for s in `cat killed15_3`; do echo sbatch $s; done > resubmit
bash resubmit

# clean up
rm -f killed15 killed152 killed15_3 resubmit 
