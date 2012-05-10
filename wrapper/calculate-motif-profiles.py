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

net = readNetwork(options.Filename)
runMotifProfileAnalysis(net,options.MotifSize,options.Bipartite,options.Weighted)
