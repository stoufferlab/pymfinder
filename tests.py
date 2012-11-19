#!/usr/bin/python

from wrapper import *

##############################################################
##############################################################
# C'est le debut
##############################################################
##############################################################

filename = "test.net"
structure = False

# calculate the network's motif structure
if structure:
    mm = motif_structure(filename,
                         motifsize = 3,
                         nrandomizations = 10,
                         usemetropolis = False,
                         )

    # print out the results
    print_motif_structure(mm,sep='\t',header=True)

# calculate the network's (node-by-node) motif participation
else:
    pp = motif_participation(filename,
                             motifsize = 3,
                             randomize = False,
                             usemetropolis = False,
                             )

    # print out the results
    print_participation(pp,sep='\t',header=True)

##############################################################
##############################################################
# C'est fini
##############################################################
##############################################################
