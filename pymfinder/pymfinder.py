#!/usr/bin/python

import mfinder.mfinder as mfinder
import sys
from roles import * #Now storing motif ID's, roles offsite.

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

# turn an mfinder-style network into a human-intelligible one
def human_network_setup(network):
    if type(network) == type("hello world"):
        # DEBUG: if we want to use a filename we need to run a check here to make sure that the node labels are integers and that there are weights
        # web.Filename = network
        network = read_links(network)
        return relabel_nodes(network)
    elif type(network) == type([1,2,3]):
        return relabel_nodes(network)
    else:
        sys.stderr.write("Whatcha talkin' 'bout Willis?\n")
        sys.exit()

# if we've relabeled the nodes, make sure the output corresponds to the input labels
def decode_stats(stats,node_dictionary):
    reverse_dictionary = dict([(j,i) for i,j in node_dictionary.items()])
    return dict([(reverse_dictionary[i],j) for i,j in stats.items()])

# if we've relabeled the nodes, make sure the output corresponds to the input labels
def decode_link_stats(stats,node_dictionary):
    reverse_dictionary = dict([(j,i) for i,j in node_dictionary.items()])
    return dict([((reverse_dictionary[i[0]],reverse_dictionary[i[1]]),j) for i,j in stats.items()])

##############################################################
##############################################################
# MOTIF GENERATING CODE
##############################################################
##############################################################

def list_motifs(motifsize,
                ):

    motifs = mfinder.list_motifs(motifsize)

    all_motifs = []
    motif_result = motifs.l
    while (motif_result != None):
        all_motifs.append(motif_result.val)
                
        motif_result = motif_result.next

    return all_motifs

def print_motifs(motifs,motifsize,outFile=None,sep=" ",links=False):
    if outFile:
        fstream = open(outFile,'w')
    else:
        fstream = sys.stdout

    for m in motifs:
        output = sep.join(["%i" % m,
                           ])

        fstream.write(output + '\n')

        if links:
            motif_edges = mfinder.motif_edges(m,motifsize)
            edge_result = motif_edges.l
            while (edge_result != None):
                edge = mfinder.get_edge(edge_result.p)
                s = int(edge.s)
                t = int(edge.t)
                output = sep.join(["%i" % s,
                                   "%i" % t,
                                   ])
                fstream.write(output + '\n')
                edge_result = edge_result.next

    if outFile:
        fstream.close()

    return   

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

def print_random_network(edges,outFile=None,sep=" ",weighted=False,header=False):
    if outFile:
        fstream = open(outFile,'w')
    else:
        fstream = sys.stdout

    if header:
        if weighted:
          output = sep.join(['target',
                             'source',
                             'weight',])
        else:
          output = sep.join(['target',
                             'source',])

        fstream.write(output + '\n')

    for trg,src,w in sorted(edges):
        if weighted:
          output = sep.join(["%i" % trg,
                             "%i" % src,
                             "%i" % w,
                             ])
        else:
          output = sep.join(["%i" % trg,
                             "%i" % src,
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
            sys.stderr.write("upping the ante bitches!\n")
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
        return participation    #participation[motif_ID][species] = number of times that species is in that motif

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
        return decode_stats(p_stats,node_dict)
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
# LINK PARTICIPATION CODE
##############################################################
##############################################################

def link_participation_stats(mfinderi,network,stoufferIDs):
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
            sys.stderr.write("upping the ante bitches!\n")
            mfinderi.MaxMembersListSz = max_count + 1
            results = mfinder.motif_participation(mfinderi)

        else:
            break

    _network = [(i,j) for i,j,k in network]

    linkparticipation = {}
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
              for n in py_members:         #n is a species
                if m!=n:
                  link=(m, n)
                  if link in _network:

                    if link not in linkparticipation:
                      linkparticipation[link] = {}
                    try:
                      linkparticipation[link][id] += 1
                    except KeyError:
                      linkparticipation[link][id] = 1
            am_l = am_l.next

        r_l = r_l.next
    r_l = results.l

    while (r_l != None):
        motif = mfinder.get_motif(r_l.p)
        id = int(motif.id)

        for link in linkparticipation:
            try:
                x = linkparticipation[link][id]
            except KeyError:
                linkparticipation[link][id] = 0

        r_l = r_l.next 

    return linkparticipation

def link_participation(network,
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

    network, node_dict = human_network_setup(network)
    pl_stats = link_participation_stats(web,network,stoufferIDs)

    try:
        return decode_link_stats(pl_stats,node_dict)
    except UnboundLocalError:
        return pl_stats

    #Dealing with (a, b), (b, a) issues:
    netfile=open(network)
    intlist=[]
    for line in netfile:
      (pred,prey)=line.split()
      pair=(int(pred),int(prey))
      pair2=(int(prey),int(pred))
      if pair2 not in intlist:
        intlist.append(pair)

    link_stats={}
    for link in pl_stats:
      if link in intlist:
        link_stats[link]={}
    for key in link_stats:
      link_stats[key]=pl_stats[key]
    return link_stats    

def print_link_participation(link_participation_stats,outFile=None,sep=" ",header=False):
    if outFile:
        fstream = open(outFile,'w')
    else:
        fstream = sys.stdout

    if header:
        output = sep.join(["link"]+list(map(str,sorted(link_participation_stats[link_participation_stats.keys()[0]].keys()))))
        fstream.write(output + '\n')

    for (a, b) in sorted(link_participation_stats.keys()):
        output = sep.join(['('+str(a)+','+str(b)+')'] + list(map(str,[j for i,j in sorted(link_participation_stats[(a, b)].items())])))
        fstream.write(output + '\n')

    if outFile:
        fstream.close()

    return

##############################################################
##############################################################
# MOTIF ROLES CODE
##############################################################
##############################################################

def role_stats(mfinderi,network,stoufferIDs,networktype):
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
            sys.stderr.write("upping the ante bitches!\n")
            mfinderi.MaxMembersListSz = max_count + 1
            results = mfinder.motif_participation(mfinderi)

        else:
            break

    _network = [(i,j) for i,j,k in network]

    possible_roles = []    
    if networktype == "unipartite":
      for motif,roles in UNIPARTITE_ROLES[mfinderi.MotifSize]:
          possible_roles += [tuple([motif] + list(role)) for role in roles]
    elif networktype == "bipartite":
      for motif,roles in BIPARTITE_ROLES[mfinderi.MotifSize]:
          possible_roles += [tuple([motif] + list(role)) for role in roles]

    roles = {}
    r_l = results.l
    members = mfinder.intArray(mfinderi.MotifSize)
    while (r_l != None):
        motif = mfinder.get_motif(r_l.p)
        id = int(motif.id)

        if networktype == "unipartite":
          mindex = [i for i,j in UNIPARTITE_ROLES[mfinderi.MotifSize]].index(id)
        elif networktype == "bipartite":
          mindex = [i for i,j in BIPARTITE_ROLES[mfinderi.MotifSize]].index(id)
        
        am_l = motif.all_members.l
        while (am_l != None):
            mfinder.get_motif_members(am_l.p, members, mfinderi.MotifSize)
            py_members = [int(members[i]) for i in xrange(mfinderi.MotifSize)]

            for m in py_members:
                if m not in roles:
                    roles[m] = {}

                npred  = sum([(othernode,m) in _network for othernode in py_members if othernode != m])
                nprey = sum([(m,othernode) in _network for othernode in py_members if othernode != m])

                key = (id, npred, nprey)
                #key = (id, indegree, outdegree)

                if key not in possible_roles:
                    if npred > 0:
                        connected_to = [othernode for othernode in py_members if othernode != m and (othernode,m) in _network]
                        npreys = [sum([(i,j) in _network for j in py_members if j != i]) for i in connected_to]
                        npreys.sort()
                        key = tuple(list(key) + [tuple(npreys)])
                    else:
                        connected_to = [othernode for othernode in py_members if othernode != m and (m,othernode) in _network]
                        npreds = [sum([(j,i) in _network for j in py_members if j != i]) for i in connected_to]
                        npreds.sort()
                        key = tuple(list(key) + [tuple(npreds)])

                if key not in possible_roles:
                    print >> sys.stderr, key
                    print >> sys.stderr, "Apparently there is a role you aren't accounting for in roles.py."
                    sys.exit()

                try:
                    roles[m][key] += 1
                except KeyError:
                    roles[m][key] = 1

            am_l = am_l.next

        r_l = r_l.next
    for n in roles:
        for r in possible_roles:
            try:
                x = roles[n][r]
            except KeyError:
                roles[n][r] = 0
    return roles

def motif_roles(network,
                motifsize = 3,
                maxmemberslistsz = 1000,
                randomize = False,
                usemetropolis = False,
                stoufferIDs = False,
                networktype="unipartite"
                ):

    print networktype
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
    
    network, node_dict = human_network_setup(network)    
    r_stats = role_stats(web,network,stoufferIDs,networktype)

    try:
        return decode_stats(r_stats,node_dict)
    except UnboundLocalError:
        return r_stats

def print_roles(role_stats,motifsize,outFile=None,sep=" ",header=False):
    if outFile:
        fstream = open(outFile,'w')
    else:
        fstream = sys.stdout

    roles = []
    for motif,rs in UNIPARTITE_ROLES[motifsize]:
        roles += [tuple([motif] + list(r)) for r in rs]

    if header:
        output = sep.join(["node"]+list(map(str,range(1,len(roles)+1))))
        fstream.write(output + '\n')

    for n in sorted(role_stats.keys()):
        output = sep.join([str(n)] + list(map(str,[role_stats[n][role] for role in roles])))
        fstream.write(output + '\n')

    if outFile:
        fstream.close()

    return

##############################################################
##############################################################
# LINK ROLES CODE
##############################################################
##############################################################

def link_role_stats(mfinderi,network,stoufferIDs): #End product should be a nested dict [link][linkrole]: count
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
            sys.stderr.write("upping the ante bitches!\n")
            mfinderi.MaxMembersListSz = max_count + 1
            results = mfinder.motif_participation(mfinderi)

        else:
            break

    _network = [(i,j) for i,j,k in network]

    possible_links = []
    for motif,links in UNIPARTITE_LINKS[mfinderi.MotifSize]:
        possible_links += [tuple([motif] + list(link)) for link in links]

    lroles = {}           #roles[link]=(motifID, npred, nprey)
    for link in _network:
      (a,b)=link
      if int(a)!=int(b):
        lroles[link]={}

    r_l = results.l
    members = mfinder.intArray(mfinderi.MotifSize)
    while (r_l != None):
        motif = mfinder.get_motif(r_l.p)
        id = int(motif.id)
        lindex = [i for i,j in UNIPARTITE_LINKS[mfinderi.MotifSize]].index(id)
        am_l = motif.all_members.l
        while (am_l != None):
            mfinder.get_motif_members(am_l.p, members, mfinderi.MotifSize)
            py_members = [int(members[i]) for i in xrange(mfinderi.MotifSize)]
            for link in lroles:  #This will avoid the (a,b) vas (b,a) difficulty.
              (sp1,sp2)=link
              if sp1 in py_members and sp2 in py_members:
                npred1  = sum([(othernode,sp1) in _network for othernode in py_members if othernode != sp1])
                nprey1 = sum([(sp1,othernode) in _network for othernode in py_members if othernode != sp1])
                npred2  = sum([(othernode,sp2) in _network for othernode in py_members if othernode != sp2])
                nprey2 = sum([(sp2,othernode) in _network for othernode in py_members if othernode != sp2])
                key = (id, (npred1, nprey1),(npred2, nprey2))

                if key not in possible_links:
                  if sp1==sp2: #Don't care about cannibalistic links
                    #print >> sys.stderr, "We don't cotton to no cannibals 'round here"
                    pass
                  else:
                    altkey = (id, (npred2, nprey2),(npred1, nprey1))
                    if altkey not in possible_links:
                      #print key, altkey, link
                      print >> sys.stderr, "Apparently there is a role you aren't accounting for in roles.py."
                      sys.exit()
                    else:
                      key=altkey
                try:
                  lroles[link][key] += 1
                except KeyError:
                  lroles[link][key] = 1

            am_l = am_l.next

        r_l = r_l.next
    for l in lroles:
        for r in possible_links:
            try:
                x = lroles[l][r]
            except KeyError:
                lroles[l][r] = 0
    return lroles

def link_roles(network,
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
    
    network, node_dict = human_network_setup(network)
    r_stats = link_role_stats(web,network,stoufferIDs)

    try:
        return decode_link_stats(r_stats,node_dict)
    except UnboundLocalError:
        return r_stats

def print_link_roles(link_role_stats,motifsize,outFile=None,sep=" ",header=False):
    if outFile:
        fstream = open(outFile,'w')
    else:
        fstream = sys.stdout

    lroles = []
    for motif,rs in UNIPARTITE_LINKS[motifsize]:
        lroles += [tuple([motif] + list(r)) for r in rs]

    if header:
        output = sep.join(["link"]+list(map(str,range(1,len(lroles)+1))))
        fstream.write(output + '\n')

    for (a, b) in sorted(link_role_stats.keys()):
        output = sep.join(['('+str(a)+','+str(b)+')'] + list(map(str,[link_role_stats[(a, b)][role] for role in lroles])))
        fstream.write(output + '\n')

    if outFile:
        fstream.close()

    return


##############################################################
##############################################################
# C'est fini
##############################################################
##############################################################