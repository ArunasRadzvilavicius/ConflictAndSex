# Calculate fixation probabilities from
# .evo files containing the  total number of 
# simulation runs and the number of fixations
# Use: ./CalcStats.py 0 -1 ListOfFilenames
#! /usr/bin/env python
import sys
import numpy
start=int(sys.argv[1])
end=int(sys.argv[2])
fnames = sys.argv[3:]
for fname in fnames:
	a = numpy.loadtxt(fname)
	mu = fname[start:end]
	ns = sum(a[:,0])
	nt = sum(a[:,1])
	print mu, 1.0*ns/nt, nt
