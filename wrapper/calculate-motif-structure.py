#!/usr/bin/python

# system libraries
import sys

# local libraries
from tools import *
from option_parser import motif_structure as option_parser

####################################
### Start of main of program
####################################

(options, args) = option_parser()

net = readNetwork(options.Filename)
runMotifAnalysis(net,options.MotifSize,options.NRandomizations)
