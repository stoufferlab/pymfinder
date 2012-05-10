#!/usr/bin/python

# system libraries
#import sys
from optparse import OptionParser

###########################################

def motif_structure():
	usage = "usage: %prog [whatever] sucker"
        parser = OptionParser(usage)

        parser.add_option("-f", "--filename", 
                          action="store",
                          dest="Filename",
                          help="name of network input file [default: stdin]",
                          )

        parser.add_option("-m", "--motifsize", 
                          action="store",
			  type="int",
                          dest="MotifSize",
                          help="calculate statistics for m-species motifs [default: 3]",
                          )

        parser.add_option("-r", "--randomizations", 
                          action="store",
			  type="int",
                          dest="NRandomizations",
                          help="compare empirical motif frequencies to an ensemble of r randomizations [default: 250]",
                          )

        parser.set_defaults(Filename = '-',
                            MotifSize = 3,
			    NRandomizations = 250,
                            )

        (options, args) = parser.parse_args()

        return (options, args)

def motif_profiles():
	usage = "usage: %prog [whatever] sucker"
        parser = OptionParser(usage)

        parser.add_option("-f", "--filename", 
                          action="store",
                          dest="Filename",
                          help="name of network input file [default: stdin]",
                          )

        parser.add_option("-m", "--motifsize", 
                          action="store",
			  type="int",
                          dest="MotifSize",
                          help="calculate statistics for m-species motifs [default: 3]",
                          )

        parser.add_option("-n", "--normalized", 
                          action="store_false",
                          dest="Normalized",
                          help="normalize the role frequencies within each motif size group [default: Yes]",
                          )

        parser.add_option("-w", "--weighted", 
                          action="store_true",
                          dest="Weighted",
                          help="calculate weighted motif profiles [default: No]",
                          )

	parser.add_option("-b", "--bipartite",
			  action="store_true",
			  dest="Bipartite",
			  help="treat the network as a bipartite network [default: No]",
			  )

        parser.set_defaults(Filename = '-',
			    MotifSize = 3,
                            Normalized = True,
                            Weighted = False,
			    Bipartite = False,
                            )

        (options, args) = parser.parse_args()

        return (options, args)

