
import mfinder.mfinder as cmfinder
import sys
from itertools import combinations

from roles import *
from datatypes import *

##############################################################
##############################################################
# GENERAL UTILITIES
##############################################################
##############################################################

def read_links(filename):
    inFile = open(filename, 'r')
    links = []
    for i in inFile.readlines():
        l = i.strip().split()
        if len(l) > 3 or len(l) < 2:
            inFile.close()
            sys.stderr.write("There is something peculiar about one of the interactions in your input file.\n")
            sys.exit()
        elif len(l) == 2:
            links += [tuple(l + [1])]
        else:
            links += [tuple(l)]

    inFile.close()

    return links

# turn any type of node label into integers (mfinder is finicky like that)
def relabel_nodes(links,stats=None):
    node_dict = {}
    for i in range(len(links)):
        try:
            s,t,w = links[i]
            w = float(w)
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
            node_dict[s] = len(node_dict)+1
        if t not in node_dict:
            node_dict[t] = len(node_dict)+1

        links[i] = (node_dict[s], node_dict[t], w)

        if stats:
            stats.add_link((s,t))
            try:
                x = stats.nodes[s]
            except KeyError:
                stats.add_node(s)
            try:
                x = stats.nodes[t]
            except KeyError:
                stats.add_node(t)

    return links, node_dict
    
def gen_mfinder_network(links):
    edges = cmfinder.intArray(len(links)*3+1)
    for i in range(len(links)):
        try:
            s,t,w = links[i]
            w = int(round(w))
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
        stats = NetworkStats()
        network = read_links(network)
        network, node_dict = relabel_nodes(network,stats)
        edges, numedges = gen_mfinder_network(network)
        return network, stats, edges, numedges, node_dict
    elif type(network) == type([1,2,3]):
        stats = NetworkStats()
        network, node_dict = relabel_nodes(network,stats)
        edges, numedges = gen_mfinder_network(network)
        return network, stats, edges, numedges, node_dict
    elif type(network) == NetworkStats:
        links, node_dict = relabel_nodes(network.links.keys())
        edges, numedges = gen_mfinder_network(links)
        return links, network, edges, numedges, node_dict
    else:
        sys.stderr.write("Uncle Sam frowns upon tax cheats.\n")
        sys.exit()

# TODO: This is a function that is not really necessary because the C code should provide the adjacency

def adjacency(links, size):
    adj=[[0.0 for x in range(size)] for y in range(size)]
    for i in links:
        adj[i[0]-1][i[1]-1]=i[2]
    return adj

# TODO: there must be better ways to write this.
def default_fweight(x):
    return sum(x)/len(x)

# if we've relabeled the nodes, make sure the output corresponds to the input labels
def decode_net(edges,node_dictionary):
    reverse_dictionary = dict([(j,i) for i,j in node_dictionary.items()])
    return [(reverse_dictionary[i],reverse_dictionary[j],k) for i,j,k in edges]

##############################################################
##############################################################
# MOTIF GENERATING CODE
##############################################################
##############################################################

def list_motifs(motifsize):

    motifs = cmfinder.list_motifs(motifsize)

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
            motif_edges = cmfinder.motif_edges(m,motifsize)
            edge_result = motif_edges.l
            while (edge_result != None):
                edge = cmfinder.get_edge(edge_result.p)
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
    web = cmfinder.mfinder_input()

    # setup the network info
    network, stats, web.Edges, web.NumEdges, node_dict = mfinder_network_setup(network)

    # parameterize the analysis
    if not usemetropolis:
        web.UseMetropolis = 0
    else:
        web.UseMetropolis = 1

    return decode_net(randomized_network(web), node_dict)
        
def randomized_network(mfinderi):
    results = cmfinder.random_network(mfinderi)
    
    edges = []
    edge_result = results.l
    while (edge_result != None):
        edge = cmfinder.get_edge(edge_result.p)
        s = int(edge.s)
        t = int(edge.t)
        w = int(edge.weight)
        edges.append((s,t,w))

        edge_result = edge_result.next

    return edges

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
                    weighted = False,
                    fweight = None
                    ):

    #if weighted and nrandomizations > 0:
    #    sys.stderr.write("Sorry, randomizations for weighted networks aren't implemented yet.\n")
    #    sys.exit()

    # initialize the heinous input struct
    web = cmfinder.mfinder_input()

    # setup the network info
    network, stats, web.Edges, web.NumEdges, node_dict = mfinder_network_setup(network)

    # add or check some basics of stats object
    if stats.motifsize:
        if stats.motifsize != motifsize:
            sys.stderr.write("You're trying to mix motif sizes.\n")
            sys.exit()
    else:
        stats.motifsize = motifsize

	# TODO: Double check this!
    if stats.networktype:
        stats.networktype = "unipartite"

    if stats.weighted!=None:
        if stats.weighted != weighted:
            sys.stderr.write("You're trying to mix two different motif analyses.\n")
            stats.weighted = weighted
    else:
        stats.weighted = weighted
        if stats.weighted:
            stats.adj = adjacency(network, len(node_dict))

    if fweight==None:
        fweight = default_fweight

    # parameterize the analysis
    web.MotifSize = motifsize
    web.NRandomizations = nrandomizations
    if not usemetropolis:
        web.UseMetropolis = 0
    else:
        web.UseMetropolis = 1

    #check if this function has already been run
    if len(stats.motifs) != 0:
        stats.motifs = dict()

    # determine all nodes' role statistics
    if stats.weighted:
        if len(stats.motifs) == 0:
            #TODO why is there a false with stoufferID
            motif_stats(web,stats,stoufferIDs=False)
        web.MaxMembersListSz = max([stats.motifs[x].real for x in stats.motifs])+1
        return weighted_motif_stats(web, network, stats, stoufferIDs, fweight)
    else:
        return motif_stats(web, stats, stoufferIDs)

        
def motif_stats(mfinderi,motif_stats,stoufferIDs):
    results = cmfinder.motif_structure(mfinderi)

    if results:
        motif_result = results.l
        while (motif_result != None):
            motif = cmfinder.get_motif_result(motif_result.p)

            motif_id = int(motif.id)
            
            motif_stats.add_motif(motif_id)
            motif_stats.motifs[motif_id].real = int(motif.real_count)
            motif_stats.motifs[motif_id].random_m = float(motif.rand_mean)
            motif_stats.motifs[motif_id].random_sd = float(motif.rand_std_dev)
            motif_stats.motifs[motif_id].real_z = float(motif.real_zscore)
            motif_stats.motifs[motif_id].weighted = 0

            motif_result = motif_result.next

    cmfinder.list64_free_mem(results)

    if stoufferIDs:
        motif_stats.use_stouffer_IDs()

    return motif_stats

def weighted_motif_stats(mfinderi, network, motif_stats, stoufferIDs, fweight):

    results = cmfinder.motif_participation(mfinderi)

    _network = set([(i,j) for i,j,k in network])

    r_l = results.l
    members = cmfinder.intArray(mfinderi.MotifSize)
    while (r_l != None):
        motif = cmfinder.get_motif(r_l.p)
        id = int(motif.id)

        am_l = motif.all_members.l
        while (am_l != None):
            cmfinder.get_motif_members(am_l.p, members, mfinderi.MotifSize)
            py_members = [int(members[i]) for i in xrange(mfinderi.MotifSize)]

            # TODO optimize next line. Ask Daniel for getting edges of a motif
            py_motif = set([(i,j) for i,j in _network if (i in py_members and j in py_members and i!=j)])

            # TODO add functional to choose the way to add the weights
            weight = [motif_stats.adj[m[0]-1][m[1]-1] for m in py_motif]

            motif_stats.motifs[id].weighted += fweight(weight)

            am_l = am_l.next

        r_l = r_l.next

    cmfinder.res_tbl_mem_free_single(results)

    if stoufferIDs:
        motif_stats.use_stouffer_IDs()

    return motif_stats


##############################################################
##############################################################
# MOTIF PARTICIPATION CODE
##############################################################
##############################################################

def motif_participation(network,
                        links = False,
                        motifsize = 3,
                        randomize = False,
                        usemetropolis = False,
                        stoufferIDs = False,
                        allmotifs = False,
                        weighted = False,
                        fweight = None
                        ):


    # do we want to randomize the network first?
    if randomize:
        #This will restart the whole object
        network = random_network(network, usemetropolis = usemetropolis)

    # initialize the heinous input struct
    web = cmfinder.mfinder_input()

    # setup the network info
    network, stats, web.Edges, web.NumEdges, node_dict = mfinder_network_setup(network)

    # add or check some basics of stats object
    if stats.motifsize:
        if stats.motifsize != motifsize:
            sys.stderr.write("You're trying to mix motif sizes.\n")
            sys.exit()
    else:
        stats.motifsize = motifsize

    if stats.networktype:
        stats.networktype = "unipartite"

    if stats.weighted!=None:
        if stats.weighted != weighted:
            sys.stderr.write("You're trying to mix two different motif analyses (weighted and not weighted). Be careful!\n")
            stats.weighted = weighted
    else:
        stats.weighted = weighted
        if stats.weighted:
            stats.adj = adjacency(network, len(node_dict))
            stats.node_dict=node_dict

    if fweight==None:
        fweight = default_fweight

    # parameterize the analysis
    web.MotifSize = motifsize
    web.Randomize = 0
    web.UseMetropolis = 0
    if len(stats.motifs) == 0:
        web.NRandomizations = 0
        web.UseMetropolis = 0
        #TODO why is there a false with stoufferID
        motif_stats(web,stats,stoufferIDs=False)

    web.MaxMembersListSz = max([stats.motifs[x].real for x in stats.motifs])+1

    #TODO should I run this here? Alternatively, I can do it inside participation
    #if stats.weighted:
    #    weighted_motif_stats(web,stats,stoufferIDs=False)

    #check if this function has already been run
    if len(stats.nodes[stats.nodes.keys()[0]].motifs) != 0:
        for x in stats.nodes.keys():
            stats.nodes[x].motifs = dict()
        if len(stats.links[stats.links.keys()[0]].motifs) != 0:
            for x in stats.links.keys():
                stats.links[x].motifs = dict()

    return participation_stats(web,network,stats,node_dict,links,stoufferIDs,allmotifs,fweight)


def participation_stats(mfinderi, network, participation, node_dict, links, stoufferIDs, allmotifs, fweight):

    results = cmfinder.motif_participation(mfinderi)

    node_dict = dict((v,k) for k,v in node_dict.iteritems())

    _network = set([(i,j) for i,j,k in network])

    possible_motifs = set(STOUFFER_MOTIF_IDS.keys())
    actual_motifs = set([])

    r_l = results.l
    members = cmfinder.intArray(mfinderi.MotifSize)
    while (r_l != None):
        motif = cmfinder.get_motif(r_l.p)
        id = int(motif.id)
        actual_motifs.add(id)

        am_l = motif.all_members.l
        while (am_l != None):
            cmfinder.get_motif_members(am_l.p, members, mfinderi.MotifSize)
            py_members = [int(members[i]) for i in xrange(mfinderi.MotifSize)]

            if participation.weighted:
                # TODO optimize next line. Ask Daniel for getting edges of a motif
                py_motif = set([(i,j) for i,j in _network if (i in py_members and j in py_members and i!=j)])

                # TODO add functional to choose the way to add the weights
                weight = [participation.adj[m[0]-1][m[1]-1] for m in py_motif]

                for m in py_members:
                    try:
                        participation.nodes[node_dict[m]].motifs[id] += 1
                    except KeyError:
                        participation.nodes[node_dict[m]].motifs[id] = 1
                    try:
                        participation.nodes[node_dict[m]].weighted_motifs[id] += fweight(weight)
                    except KeyError:
                        participation.nodes[node_dict[m]].weighted_motifs[id] = fweight(weight)
            else:
                for m in py_members:
                    try:
                        participation.nodes[node_dict[m]].motifs[id] += 1
                    except KeyError:
                        participation.nodes[node_dict[m]].motifs[id] = 1

            if links:
                for m, n in combinations(py_members, 2):
                    if (node_dict[m], node_dict[n]) in participation.links:
                        try:
                            participation.links[(node_dict[m], node_dict[n])].motifs[id] += 1
                        except KeyError:
                            participation.links[(node_dict[m], node_dict[n])].motifs[id] = 1

                        if participation.weighted:
                            try:
                                participation.links[(node_dict[m], node_dict[n])].weighted_motifs[id] += participation.adj[m-1][n-1]
                            except KeyError:
                                participation.links[(node_dict[m], node_dict[n])].weighted_motifs[id] = participation.adj[m-1][n-1]


                    if (node_dict[n], node_dict[m]) in participation.links:
                        try:
                            participation.links[(node_dict[n], node_dict[m])].motifs[id] += 1
                        except KeyError:
                            participation.links[(node_dict[n], node_dict[m])].motifs[id] = 1

                        if participation.weighted:
                            try:
                                participation.links[(node_dict[n], node_dict[m])].weighted_motifs[id] += participation.adj[n-1][m-1]
                            except KeyError:
                                participation.links[(node_dict[n], node_dict[m])].weighted_motifs[id] = participation.adj[n-1][m-1]


            am_l = am_l.next

        r_l = r_l.next

    cmfinder.res_tbl_mem_free_single(results)


    if not allmotifs:
        possible_motifs = actual_motifs


    for r in possible_motifs:
        for n in participation.nodes:
            try:
                x = participation.nodes[n].motifs[r]
            except KeyError:
                participation.nodes[n].motifs[r] = 0

            if participation.weighted:
                try:
                    x = participation.nodes[n].weighted_motifs[r]
                except KeyError:
                    participation.nodes[n].weighted_motifs[r] = 0

        if links:
            for n in participation.links:
                try:
                    x = participation.links[n].motifs[r]
                except KeyError:
                    participation.links[n].motifs[r] = 0

            if participation.weighted:
                try:
                    x = participation.links[n].weighted_motifs[r]
                except KeyError:
                    participation.links[n].weighted_motifs[r] = 0


    if stoufferIDs:
        participation.use_stouffer_IDs()
        
    return participation


##############################################################
##############################################################
# MOTIF ROLES CODE
##############################################################
##############################################################

def motif_roles(network,
                links=False,
                motifsize = 3,
                randomize = False,
                usemetropolis = False,
                stoufferIDs = False,
                networktype = "unipartite",
                allroles = False
                ):


    # do we want to randomize the network first?
    if randomize:
        #This will restart the whole object
        network = random_network(network, usemetropolis = usemetropolis)

    # initialize the heinous input struct
    web = cmfinder.mfinder_input()

    # setup the network info
    network, stats, web.Edges, web.NumEdges, node_dict = mfinder_network_setup(network)

    # add or check some basics of stats object
    if stats.motifsize:
        if stats.motifsize != motifsize:
            sys.stderr.write("You're trying to mix motif sizes.\n")
            sys.exit()
    else:
        stats.motifsize = motifsize

    if stats.networktype:
        if stats.networktype != networktype:
            sys.stderr.write("You're trying to mix network types.\n")
            sys.exit()
    else:
        stats.networktype = networktype

    # parameterize the analysis
    web.MotifSize = motifsize
    web.Randomize = 0
    web.UseMetropolis = 0
    if len(stats.motifs) == 0:
        web.NRandomizations = 0
        web.UseMetropolis = 0
        motif_stats(web,stats,stoufferIDs=False)

    web.MaxMembersListSz = max([stats.motifs[x].real for x in stats.motifs])+1

    #check if this function has already been run
    if len(stats.nodes[stats.nodes.keys()[0]].roles) != 0:
        for x in stats.nodes.keys():
            stats.nodes[x].roles = dict()
        if len(stats.links[stats.links.keys()[0]].roles) != 0:
            for x in stats.links.keys():
                stats.links[x].roles = dict()

    # determine all nodes' role statistics
    return role_stats(web,stats,node_dict,network,links,networktype,stoufferIDs,allroles)


def role_stats(mfinderi,roles,node_dict,network,links,networktype,stoufferIDs,allroles):

    results = cmfinder.motif_participation(mfinderi)

    node_dict = dict((v,k) for k,v in node_dict.iteritems())

    possible_roles = set([])
    actual_roles = set([])
    if networktype == "unipartite":
      for m,r in UNIPARTITE_ROLES[mfinderi.MotifSize]:
          possible_roles.update([tuple([m] + list(x)) for x in r])
    elif networktype == "bipartite":
      for m,r in BIPARTITE_ROLES[mfinderi.MotifSize]:
          possible_roles.update([tuple([m] + list(x)) for x in r])

    if links:
        possible_linkroles = set([])
        actual_linkroles = set([])
        for m,l in UNIPARTITE_LINKS_ROLES[mfinderi.MotifSize]:
            possible_linkroles.update([tuple([m] + list(x)) for x in l])


    _network = set([(i,j) for i,j,k in network])
    
    r_l = results.l
    members = cmfinder.intArray(mfinderi.MotifSize)
    while (r_l != None):
        motif = cmfinder.get_motif(r_l.p)
        id = int(motif.id)
        
        am_l = motif.all_members.l
        while (am_l != None):
            cmfinder.get_motif_members(am_l.p, members, mfinderi.MotifSize)
            py_members = [int(members[i]) for i in xrange(mfinderi.MotifSize)]

            #TODO remove self-edges and be coherent below
            py_motif = set([(i,j) for i,j in _network if (i in py_members and j in py_members)])

            for m in py_members:

                npred  = sum([(othernode,m) in py_motif for othernode in py_members if othernode != m])
                nprey = sum([(m,othernode) in py_motif for othernode in py_members if othernode != m])

                key = (id, npred, nprey)

                # if the node's in and out degrees are insufficient to discern its role
                # we will add the degrees of the nodes it interacts with (its neighbors)
                if key not in possible_roles:
                    if npred > 0:
                        connected_to = set([othernode for othernode in py_members if othernode != m and (othernode,m) in py_motif])
                        npreys = [sum([(i,j) in py_motif for j in py_members if j != i]) for i in connected_to]
                        npreys.sort()
                        key = tuple(list(key) + [tuple(npreys)])
                    else:
                        connected_to = set([othernode for othernode in py_members if othernode != m and (m,othernode) in py_motif])
                        npreds = [sum([(j,i) in py_motif for j in py_members if j != i]) for i in connected_to]
                        npreds.sort()
                        key = tuple(list(key) + [tuple(npreds)])

                if key not in possible_roles:
                    print >> sys.stderr, key
                    print >> sys.stderr, "Apparently there is a role you aren't accounting for in roles.py."
                    sys.exit()

                try:
                    roles.nodes[node_dict[m]].roles[key] += 1
                except KeyError:
                    roles.nodes[node_dict[m]].roles[key] = 1

                actual_roles.add(key)

                if links:
                    for n in py_members:
                        if n == m:
                            continue

                        if (node_dict[n],node_dict[m]) in roles.links:
                            npred1  = sum([(othernode,n) in _network for othernode in py_members if othernode != n])
                            nprey1 = sum([(n,othernode) in _network for othernode in py_members if othernode != n])
                            npred2  = sum([(othernode,m) in _network for othernode in py_members if othernode != m])
                            nprey2 = sum([(m,othernode) in _network for othernode in py_members if othernode != m])
                            key = (id, (npred1, nprey1),(npred2, nprey2))

                            if key not in possible_linkroles:
                                key = (id, (npred2, nprey2),(npred1, nprey1))

                            if key not in possible_linkroles:
                                print >> sys.stderr, key
                                print >> sys.stderr, "Apparently there is a role you aren't accounting for in roles.py."

                            try:
                                roles.links[(node_dict[n],node_dict[m])].roles[key] += 1
                            except KeyError:
                                roles.links[(node_dict[n],node_dict[m])].roles[key] = 1

                            actual_linkroles.add(key)

            am_l = am_l.next

        r_l = r_l.next

    cmfinder.res_tbl_mem_free_single(results)

    if not allroles:
        possible_roles = actual_roles
	if links:
            possible_linkroles = actual_linkroles


    for n in roles.nodes:
        for r in possible_roles:
            try:
                x = roles.nodes[n].roles[r]
            except KeyError:
                roles.nodes[n].roles[r] = 0

    if links:
        for n in roles.links:
            for r in possible_linkroles:
                try:
                    x = roles.links[n].roles[r]
                except KeyError:
                    roles.links[n].roles[r] = 0


    if stoufferIDs:
        roles.use_stouffer_IDs()

    return roles


##############################################################
##############################################################
# RUN IT ALL!
##############################################################
##############################################################

def pymfinder(network,
              links=False,
              motifsize = 3,
              stoufferIDs = None,
              allmotifs = False,
              nrandomizations = 0,
              randomize = False,
              usemetropolis = False,
              networktype = "unipartite"
              ):


    stats = motif_participation(network,
                                links = links,
                                motifsize = motifsize,
                                randomize = randomize,
                                usemetropolis = usemetropolis,
                                stoufferIDs = stoufferIDs,
                                allmotifs = allmotifs)

    stats = motif_roles(stats,
                        links = links,
                        motifsize = motifsize,
                        randomize = False,
                        usemetropolis = usemetropolis,
                        stoufferIDs = stoufferIDs,
                        networktype = networktype,
                        allroles = allmotifs)

    if nrandomizations != 0:
        stats = motif_structure(stats,
                                motifsize = motifsize,
                                nrandomizations = nrandomizations,
                                usemetropolis = usemetropolis,
                                stoufferIDs = stoufferIDs)



    return stats

##############################################################
##############################################################
# C'est fini
##############################################################
##############################################################
