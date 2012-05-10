#!/usr/bin/python

# system libraries
import sys

# local libraries
from tools import *
from option_parser import motif_profiles as option_parser

####################################
### Start of main of program
####################################

(options, args) = option_parser()

if options.Filename == "-":
  net = readNetwork(sys.stdin)
else:
  net = readNetwork(options.Filename)

runMotifProfileAnalysis(net,options.MotifSize,options.Bipartite,options.Weighted)
