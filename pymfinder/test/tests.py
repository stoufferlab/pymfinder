#!/usr/bin/python

from pymfinder import *

##############################################################
##############################################################
# C'est le debut
##############################################################
##############################################################

import os
filename = os.path.dirname(__file__) + "./../data/test.net"

motifs = True
randomize = True
structure = True
participation = True
roles = True

# list all motifs of a specific size
if motifs:
    msize = 4
    m = list_motifs(motifsize = msize,
                    )
    print_motifs(m,
                 motifsize = msize,
                 outFile=None,
                 sep=' ',
                 links=True,
                 )

# generate a randomized version of the network
if randomize:
    rnd = random_network(filename,
                         usemetropolis = False,
                         )

    # print out the results
    print_random_network(rnd,
                         outFile=None,
                         sep=' ',
                         header=True,
                         )

# calculate the network's motif structure
if structure:
    mm = motif_structure(filename,
                         motifsize = 3,
                         nrandomizations = 250,
                         usemetropolis = False,
                         stoufferIDs = False,
                         )

    # print out the results
    print_motif_structure(mm,
                          outFile=None,
                          sep='\t',
                          header=True,
                          )

# calculate the network's (node-by-node) motif participation
if participation:
    pp = motif_participation(filename,
                             motifsize = 3,
                             randomize = False,
                             usemetropolis = False,
                             stoufferIDs = False,
                             )

    # print out the results
    print_participation(pp,
                        outFile=None,
                        sep='\t',
                        header=True)

# calculate the network's (node-by-node) motif-role participation
if roles:
    rr = motif_roles(filename,
                     motifsize = 3,
                     randomize = False,
                     usemetropolis = False,
                     stoufferIDs = False,
                     )

    # print out the results
    print_roles(rr,
                motifsize = 3,
                outFile=None,
                sep=' ',
                header=True)
##############################################################
##############################################################
# C'est fini
##############################################################
##############################################################
