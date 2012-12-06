#!/usr/bin/python

import mfinder.mfinder as mfinder
import sys

##############################################################
##############################################################
# USEFUL GLOBAL VARIABLES
##############################################################
##############################################################

STOUFFER_MOTIF_IDS = {12:  'S1',
                      38:  'S2',
                      98:  'S3',
                      36:  'S4',
                      6:   'S5',
                      46:  'D1',
                      108: 'D2',
                      14:  'D3',
                      74:  'D4',
                      102: 'D5',
                      238: 'D6',
                      110: 'D7',
                      78:  'D8',
                      }

##############################################################
##############################################################
# GENERAL UTILITIES
##############################################################
##############################################################

def read_links(filename):
    inFile = open(filename, 'r')
    links = [i.strip().split() for i in inFile.readlines()]
    inFile.close()

    if links:
        for i in range(len(links)):
            if len(links[i]) > 3 or len(links[i]) < 2:
                sys.stderr.write("There is something peculiar about one of the interactions in your input file.\n")
                sys.exit() 
            elif len(links[i]) == 2:
                links[i] = links[i] + [1]

    return [tuple(i) for i in links]

# turn any type of node label into integers (mfinder is finicky like that)
def relabel_nodes(links):
    node_dict = {}
    for i in range(len(links)):
        try:
            s,t,w = links[i]
            w = int(w)
        except ValueError:
            s,t = links[i]
            w = 1

        try:
            s = int(s)
        except:
            pass

        try:
            t = int(t)
        except:
            pass

        if s not in node_dict:
            node_dict[s] = len(node_dict) + 1
        if t not in node_dict:
            node_dict[t] = len(node_dict) + 1

        links[i] = (node_dict[s], node_dict[t], w)

    return links, node_dict
    
def gen_mfinder_network(links):
    edges = mfinder.intArray(len(links)*3+1)
    for i in range(len(links)):
        try:
            s,t,w = links[i]
            w = int(w)
        except ValueError:
            s,t = links[i]
            w = 1

        edges[3*i+1] = s
        edges[3*i+2] = t
        edges[3*i+3] = w

    return edges, len(links)

# populate the network info
def mfinder_network_setup(network):
    if type(network) == type("hello world"):
        # DEBUG: if we want to use a filename we need to run a check here to make sure that the node labels are integers and that there are weights
        # web.Filename = network
        network = read_links(network)
        network, node_dict = relabel_nodes(network)
        edges, numedges = gen_mfinder_network(network)
        return edges, numedges, node_dict
    elif type(network) == type([1,2,3]):
        network, node_dict = relabel_nodes(network)
        edges, numedges = gen_mfinder_network(network)
        return edges, numedges, node_dict
    else:
        sys.stderr.write("Uncle Sam frowns upon tax cheats.\n")
        sys.exit()

##############################################################
##############################################################
# RANDOM NETWORK CODE
##############################################################
##############################################################

def random_network(network,
                   usemetropolis = False,
                   ):

    # initialize the heinous input struct
    web = mfinder.mfinder_input()

    # setup the network info
    web.Edges, web.NumEdges, node_dict = mfinder_network_setup(network)

    # parameterize the analysis
    if not usemetropolis:
        web.UseMetropolis = 0
    else:
        web.UseMetropolis = 1

    return randomized_network(web)
        
def randomized_network(mfinderi):
    results = mfinder.random_network(mfinderi)
    
    edges = []
    edge_result = results.l
    while (edge_result != None):
        edge = mfinder.get_edge(edge_result.p)
        s = int(edge.s)
        t = int(edge.t)
        w = int(edge.weight)
        edges.append((s,t,w))

        edge_result = edge_result.next

    return edges

def print_random_network(edges,outFile=None,sep=" ",header=False):
    if outFile:
        fstream = open(outFile,'w')
    else:
        fstream = sys.stdout

    if header:
        output = sep.join(['target',
                           'source',
                           'weight',])

        fstream.write(output + '\n')

    for trg,src,w in sorted(edges):
        output = sep.join(["%i" % trg,
                           "%i" % src,
                           "%i" % w,
                           ])

        fstream.write(output + '\n')

    if outFile:
        fstream.close()

    return

##############################################################
##############################################################
# MOTIF STRUCTURE CODE
##############################################################
##############################################################

def motif_structure(network,
                    motifsize = 3,
                    nrandomizations = 0,
                    usemetropolis = False,
                    stoufferIDs = None,
                    ):

    # initialize the heinous input struct
    web = mfinder.mfinder_input()

    # setup the network info
    web.Edges, web.NumEdges, node_dict = mfinder_network_setup(network)

    # parameterize the analysis
    web.MotifSize = motifsize
    web.NRandomizations = nrandomizations
    if not usemetropolis:
        web.UseMetropolis = 0
    else:
        web.UseMetropolis = 1

    return motif_stats(web,stoufferIDs)
        
def motif_stats(mfinderi,stoufferIDs):
    results = mfinder.motif_structure(mfinderi)

    motif_stats = {}
    motif_result = results.l
    while (motif_result != None):
        motif = mfinder.get_motif_result(motif_result.p)

        motif_id = int(motif.id)
        if motif_id in motif_stats:
            sys.stderr.write("A motif has appeared twice. How odd.\n")
        else:
            motif_stats[motif_id] = dict()
            motif_stats[motif_id]['real'] = int(motif.real_count)
            motif_stats[motif_id]['rand'] = float(motif.rand_mean)
            motif_stats[motif_id]['srand'] = float(motif.rand_std_dev)
            motif_stats[motif_id]['zscore'] = float(motif.real_zscore)

        motif_result = motif_result.next

    if stoufferIDs and mfinderi.MotifSize == 3:
        return dict([(STOUFFER_MOTIF_IDS[id],motif_stats[id]) for id in motif_stats])
    else:
        return motif_stats

    return motif_stats

def print_motif_structure(motif_stats,outFile=None,sep=" ",header=False):
    if outFile:
        fstream = open(outFile,'w')
    else:
        fstream = sys.stdout

    if header:
        output = sep.join(['motif',
                           'real',
                           'rand',
                           'srand',
                           'zscore',])
        fstream.write(output + '\n')

    for m in sorted(motif_stats.keys()):
        output = sep.join(["%s" % str(m),
                           "%i" % motif_stats[m]['real'],
                           "%.3f" % motif_stats[m]['rand'],
                           "%.3f" % motif_stats[m]['srand'],
                           "%.3f" % motif_stats[m]['zscore'],
                           ])
        fstream.write(output + '\n')

    if outFile:
        fstream.close()

    return

##############################################################
##############################################################
# MOTIF PARTICIPATION CODE
##############################################################
##############################################################

def participation_stats(mfinderi,stoufferIDs):
    results = mfinder.motif_participation(mfinderi)

    maxed_out_member_list = False
    max_count = 0
    while True:
        maxed_out_member_list = False

        r_l = results.l
        while (r_l != None):
            motif = mfinder.get_motif(r_l.p)

            if(int(motif.count) != motif.all_members.size):
                maxed_out_member_list = True
                max_count = max(max_count, int(motif.count))

            r_l = r_l.next

        if maxed_out_member_list:
            sys.stderr.write("upping the anty bitches!\n")
            mfinderi.MaxMembersListSz = max_count + 1
            results = mfinder.motif_participation(mfinderi)

        else:
            break

    participation = {}
    r_l = results.l
    members = mfinder.intArray(mfinderi.MotifSize)
    while (r_l != None):
        motif = mfinder.get_motif(r_l.p)
        id = int(motif.id)

        am_l = motif.all_members.l
        while (am_l != None):
            mfinder.get_motif_members(am_l.p, members, mfinderi.MotifSize)
            py_members = [int(members[i]) for i in xrange(mfinderi.MotifSize)]

            for m in py_members:
                if m not in participation:
                    participation[m] = {}

                try:
                    participation[m][id] += 1
                except KeyError:
                    participation[m][id] = 1

            am_l = am_l.next

        r_l = r_l.next

    r_l = results.l
    while (r_l != None):
        motif = mfinder.get_motif(r_l.p)
        id = int(motif.id)

        for n in participation:
            try:
                x = participation[n][id]
            except KeyError:
                participation[n][id] = 0

        r_l = r_l.next

    if stoufferIDs and mfinderi.MotifSize == 3:
        for n in participation:
            participation[n] = dict([(STOUFFER_MOTIF_IDS[id],participation[n][id]) for id in participation[n]])
        return participation
    else:
        return participation    

def decode_participation(participation_stats,node_dictionary):
    reverse_dictionary = dict([(j,i) for i,j in node_dictionary.items()])
    return dict([(reverse_dictionary[i],j) for i,j in participation_stats.items()])

def motif_participation(network,
                        motifsize = 3,
                        maxmemberslistsz = 1000,
                        randomize = False,
                        usemetropolis = False,
                        stoufferIDs = False,
                        ):

    # initialize the heinous input struct
    web = mfinder.mfinder_input()

    # setup the network info
    web.Edges, web.NumEdges, node_dict = mfinder_network_setup(network)

    # parameterize the analysis
    web.MotifSize = motifsize
    web.MaxMembersListSz = maxmemberslistsz

    # do we want to randomize the network first?
    if not randomize:
        web.Randomize = 0
    else:
        web.Randomize = 1
        # if so, should we use the metropolis algorithm?
        if not usemetropolis:
            web.UseMetropolis = 0
        else:
            web.UseMetropolis = 1
        
    p_stats = participation_stats(web,stoufferIDs)

    try:
        return decode_participation(p_stats,node_dict)
    except UnboundLocalError:
        return p_stats

def print_participation(participation_stats,outFile=None,sep=" ",header=False):
    if outFile:
        fstream = open(outFile,'w')
    else:
        fstream = sys.stdout

    if header:
        output = sep.join(["node"]+list(map(str,sorted(participation_stats[participation_stats.keys()[0]].keys()))))
        fstream.write(output + '\n')

    for n in sorted(participation_stats.keys()):
        output = sep.join([str(n)] + list(map(str,[j for i,j in sorted(participation_stats[n].items())])))
        fstream.write(output + '\n')

    if outFile:
        fstream.close()

    return

##############################################################
##############################################################
# C'est fini
##############################################################
##############################################################
