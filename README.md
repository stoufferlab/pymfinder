pymfinder
=========

Pymfinder is a library of functions dealing with motif structure for networks and roles for species.
A wrapper python script is required to manage input files and feed them to pymfinder functions.

Input networks must be a list of tuples (pred, prey, weight).
If there are no interaction strengths available, set all weights=1.

Minimal input for motif_roles is an appropriate tuple-form network. Motifsize, etc. are set by default.
Minimal input for print_roles in output from motif_roles.
In the wrapper function calling either of these, specify whether stoufferIDs=True (3-motifs only) or leave as default (False).

Minimal input for link_roles is an appropriate tuple-form network, as motif_roles.
Minimal input for print_link_roles is output from link_roles, motifsize. Not sure why this doesn't follow the default.

