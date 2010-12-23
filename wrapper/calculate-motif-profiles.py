#!/usr/bin/python

# system libraries
import sys

# local libraries
from tools import *

####################################
### Start of main of program
####################################

filename=sys.argv[1]

if filename == '-':
    filename = sys.stdin

try:
    motifsize = int(sys.argv[2])
except:
    motifsize = 3

try:
    nrandom = int(sys.argv[3])
except:
    nrandom = 250

net = readNetwork(filename)

runMotifProfileAnalysis(net,motifsize,nrandom)

sys.exit()
