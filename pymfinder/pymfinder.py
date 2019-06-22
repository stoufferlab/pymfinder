
import mfinder.mfinder as cmfinder
import sys
from itertools import combinations, permutations, product, chain
from math import sqrt
from roles import *
from datatypes import *
import numpy as np

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
            sys.stderr.write("Error: there is something peculiar about one of the interactions in your input file.\n")
            sys.exit()
        elif len(l) == 2:
            links += [tuple(l + [1])]
        else:
            links += [tuple(l)]

    inFile.close()

    return links

def invalidlayer(s):
    try: 
        int(s)
        if int(s)>0:
            return False
        else:
            return True
    except ValueError:
        return True

def read_layers(filename, stats):
    nodes_dict = dict([(stats.nodes[i].id, i) for i in stats.nodes.keys()])
    inFile = open(filename, 'r')
    layers = dict()
    for i in inFile.readlines():
        l = i.strip().split()
        if len(l) > 2:
            inFile.close()
            sys.stderr.write("Error: there is something peculiar about one of the species in your layers file.\n")
            sys.exit()
        if invalidlayer(l[1]):
            inFile.close()
            sys.stderr.write("Error: the labelling of the layers need to consecutive integers > 0.\n")
            sys.exit()
        if l[0] in layers:
            inFile.close()
            sys.stderr.write("Error: species defined multiple times in your layers file.\n")
            sys.exit()
        if l[0] not in nodes_dict:
            inFile.close()
            sys.stderr.write("Error: a species defined in your layers file is not present in the network.\n")
            sys.exit()
        layers[l[0]]=int(l[1])

    inFile.close()
    nlayers = np.asarray(layers.values())
    if np.max(nlayers) != len(np.unique(nlayers)):
        inFile.close()
        sys.stderr.write("Error: the labelling of the layers need to consecutive integers > 0.\n")
        sys.exit()
    else:
        nlayers = np.max(nlayers)

    for i in nodes_dict:
        if i not in layers:
            inFile.close()
            sys.stderr.write("Error: a species defined in your network is not present in your layers file.\n")
            sys.exit()
        stats.nodes[nodes_dict[i]].layer=layers[i]

    return nlayers, stats


# turn any type of node label into integers (mfinder is finicky like that)
def relabel_nodes(links,stats, buildon=False):
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

        if buildon:
            links[i] = (s, t, stats.links[(s,t)].weight)

        else:
            links[i] = (node_dict[s], node_dict[t], w)
            stats.add_link(link_id = (node_dict[s], node_dict[t]), link_name = (s,t))
            stats.links[(node_dict[s], node_dict[t])].weight = w
            try:
                x = stats.nodes[node_dict[s]]
            except KeyError:
                stats.add_node(node_id = node_dict[s], node_name = s)
            try:
                x = stats.nodes[node_dict[t]]
            except KeyError:
                stats.add_node(node_id = node_dict[t], node_name = t)

    return links
    
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
        network = relabel_nodes(network,stats, buildon=False)
        edges, numedges = gen_mfinder_network(network)
        return network, stats, edges, numedges
    elif type(network) == type([1,2,3]):
        stats = NetworkStats()
        network = relabel_nodes(network,stats, buildon=False)
        edges, numedges = gen_mfinder_network(network)
        return network, stats, edges, numedges
    elif type(network) == NetworkStats:
        links = relabel_nodes(network.links.keys(), network, buildon=True)
        edges, numedges = gen_mfinder_network(links)
        return links, network, edges, numedges
    else:
        sys.stderr.write("Error: this is an invalid nework input.\n")
        sys.exit()


def default_fweight(x):
    return sum(x)/len(x)


def confidence_interval(data, confidence=0.75):
    av=np.mean(data)
    m=np.median(data)
    sd=np.std(data)
    n=len(data)
    if n==1:
        return data[0], 0.0, data[0],data[0],data[0]
    n_data=np.sort(data)
    mi=n_data[int(round(n*(1-confidence)))]
    ma=n_data[int(round(n*confidence)-1)]
    return av, sd, m, ma, mi

def check_bipartite(motif):
    bipartite = [(np.sum(motif[i,:])!=0) and (np.sum(motif[:,i])!=0) for i in range(0,len(motif))]
    return any(bipartite)

def build_motif_from_id(m, motifsize):
    motif = map(int, format(m,"b").zfill(motifsize**2))
    motif.reverse()
    motif = np.asarray(motif).reshape(motifsize,motifsize)
    return motif

def generate_key(motif, nlayers, layers_method="complete"):
    lperm = list(product(range(1,nlayers+1), repeat=len(motif)))
    nprey = np.sum(motif, 1)
    npred = np.sum(motif, 0)
    k = nprey+npred
    roles = dict()
    roles_ref = dict()
    if nlayers==1:
        for i in range(0, len(motif)):
            key = (npred[i],nprey[i])
            extra = (tuple(np.sort(nprey[motif[:,i]==1])), tuple(np.sort(npred[motif[i,:]==1])))
            try:
                x=roles[key]
            except KeyError:
                roles[key] = []
            roles[key] += [extra]

            ###Useful reference dictionary to locate roles and visualize motifs
            try:
                x=roles_ref[extra]
            except KeyError:
                roles_ref[extra] = dict()
            roles_ref[extra][i] = 1
            ###

    else:
        for j in lperm:
            _key=dict()
            _extra=dict()
            for i in range(0, len(motif)):
                combi = np.asarray(list(j))
                combi_1 = combi[motif[:,i]==1]
                combi_2 = combi[motif[i,:]==1]
                k_1 = nprey[motif[:,i]==1]
                k_2 = npred[motif[i,:]==1]
                if layers_method=="simple":
                    key = tuple([j[i]] + list((np.sum(combi_1==j[i]), np.sum(combi_2==j[i]), np.sum(combi_1!=j[i]), np.sum(combi_2!=j[i]))))
                    extra = tuple([j[i]] + list((tuple(np.sort(k_1[combi_1==j[i]])), tuple(np.sort(k_2[combi_2==j[i]])), tuple(np.sort(k_1[combi_1!=j[i]])), tuple(np.sort(k_2[combi_2!=j[i]])))))
                else:
                    key = tuple([j[i]] + list(chain.from_iterable((np.sum(combi_1==l), np.sum(combi_2==l)) for l in range(1,nlayers+1))))
                    extra = tuple([j[i]] + list(chain.from_iterable((tuple(np.sort(k_1[combi_1==l])), tuple(np.sort(k_2[combi_2==l]))) for l in range(1,nlayers+1))))

                _key[i]=key
                _extra[i]=extra

            lkey = int("".join([str(i[0]) for i in sorted(_extra.values())]))

            for i in range(0, len(motif)):
                key=tuple(list(_key[i])+[lkey])
                extra=tuple(list(_extra[i])+[lkey])

                try:
                    x=roles[key]
                except KeyError:
                    roles[key] = []
                roles[key] += [extra]

                ###Useful reference dictionary to locate roles and visualize motifs
                try:
                    x=roles_ref[extra]
                except KeyError:
                    roles_ref[extra] = dict()
                roles_ref[extra][i] = j
                ###

    _roles=[]
    for i in roles.keys():
        if len(set(roles[i]))==1:
            _roles += [i]
            ###Useful reference dictionary to locate roles and visualize motifs
            roles_ref[i] = roles_ref[list(set(roles[i]))[0]]
            ###

        else:
            for j in set(roles[i]):
                if j in roles:
                    sys.stderr.write("Check with the developer because something is wrong!\n")
                    sys.exit()
                else:
                    _roles += [j]

    return _roles, roles_ref

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

def print_motifs(motifsize,motifID=None,outFile=None,links=False,sep=" "):
    if outFile:
        fstream = open(outFile,'w')
    else:
        fstream = sys.stdout

    if motifID:
        motifs=[x for x in list_motifs(motifsize) if int(x)==motifID]
    else:
        motifs=list_motifs(motifsize)

    if motifs==[]:
        sys.stderr.write("Error: this motif does not exist.\n")
        sys.exit()

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

def print_role(motifsize, role):
    motifid = role[0]
    role = tuple(role[1:])
    nlayers = int((len(role)-1)*0.5)
    if nlayers==0:
        nlayers=1
    motif = build_motif_from_id(motifid, motifsize)
    extra , reference = generate_key(motif, nlayers)

    if role not in reference:
        sys.stderr.write("Error: this role does not exist.\n")
        sys.exit()

    reference = reference[role]

    output = ""
    for m in reference:
        if reference[m]==1:
            output += "  "+"_"*motifsize+"\n"
        else:
            output += "   "+"_"*motifsize+"\n"

        for r in range(motifsize):
            if r==m:
                output += "*"
            else:
                output += " "
            if reference[m]!=1:
                output += str(reference[m][r])
            output += "|"
            for c in range(motifsize):
                output += str(motif[r,c])
            output += "\n"
        output += "\n"

    print output

def generate_role_files(motifsize, networktype="unipartite", nlayers=1, layers_method="complete"):

    if motifsize < 2:
        sys.stderr.write("Error: this is not a valid motif size.\n")
        sys.exit()

    if motifsize > 4 and networktype=="unipartite":
        sys.stderr.write("Error: we can generate the unipartite roles for this motif size.\n")
        sys.exit()

    if motifsize > 6:
        sys.stderr.write("Error: we can generate the bipartite roles for this motif size.\n")
        sys.exit()

    #TODO: Sort out size 6
    if motifsize == 6:
        all_motifs = [545392672, 62, 4261936, 820, 4260912, 316, 4262960, 828, 4328496, 1836, 4394032, 1852,  12782640,  3900,
                   34352, 33848, 34344, 33336, 99880, 34360, 35896, 100912, 36408, 99896, 100920, 101944, 233016]
    else:
        motifs = cmfinder.list_motifs(motifsize)
        all_motifs = []
        motif_result = motifs.l
        while (motif_result != None):
            all_motifs.append(motif_result.val)
            motif_result = motif_result.next

    roles = []
    for m in all_motifs:
        motif = build_motif_from_id(m, motifsize)
        if networktype!="unipartite" and check_bipartite(motif):
            continue
        _roles , extra = generate_key(motif, nlayers, layers_method)
        roles += [(m, _roles)]

    finalform=set([])
    for m,r in roles:
        finalform.update([tuple([m] + list(x)) for x in r])

    return finalform

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
    network, stats, web.Edges, web.NumEdges = mfinder_network_setup(network)

    # parameterize the analysis
    if not usemetropolis:
        web.UseMetropolis = 0
    else:
        web.UseMetropolis = 1

    return randomized_network(web)
        
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
                    allmotifs = False,
                    stoufferIDs = False,
                    weighted = False,
                    fweight = None
                    ):

    if motifsize < 2:
        sys.stderr.write("Error: this is not a valid motif size.\n")
        sys.exit()

    if motifsize > 8:
        sys.stderr.write("Error: this is not a recommended motif size.\n")
        sys.exit()

    if motifsize > 4 and allmotifs:
        sys.stderr.write("Warning: 'allmotifs' will be ignored for this motif size and motif_structure will only register existing motifs in the real network.\n")
        allmotifs=False

    if weighted and nrandomizations > 0:
        sys.stderr.write("Warning: the analysis of weighted motifs won't be performed for the randomized networks, only for the real one. There are different ways to randomize weighted networks, you could define your own and run motif_structure multiple times to find the random distribution of weighted motifs.\n")

    # initialize the heinous input struct
    web = cmfinder.mfinder_input()

    # setup the network info
    network, stats, web.Edges, web.NumEdges = mfinder_network_setup(network)

    # add or check some basics of stats object
    if stoufferIDs and motifsize!=3:
        sys.stderr.write("Warning: 'stoufferIDs' can only be true when 'motifsize=3' in unipartite networks.\n")
        stats.stoufferIDs = False
    else:
        stats.stoufferIDs = stoufferIDs

    if stats.motifsize:
        if stats.motifsize != motifsize:
            sys.stderr.write("Error: you're trying to mix motif sizes.\n")
            sys.exit()
    else:
        stats.motifsize = motifsize

    if stats.networktype:
        stats.networktype = "unipartite"

    if stats.weighted!=None:
        if stats.weighted != weighted:
            sys.stderr.write("Error: you're trying to mix two different motif analyses.\n")

    stats.weighted = weighted

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
            motif_stats(web,stats, allmotifs)
        web.MaxMembersListSz = max([stats.motifs[x].real for x in stats.motifs])+1
        return weighted_motif_stats(web, stats, fweight, allmotifs)
    else:
        return motif_stats(web, stats, allmotifs)

        
def motif_stats(mfinderi,motif_stats, allmotifs):
    results = cmfinder.motif_structure(mfinderi)

    if results:
        motif_result = results.l
        while (motif_result != None):
            motif = cmfinder.get_motif_result(motif_result.p)

            motif_id = int(motif.id)
            
            if allmotifs or int(motif.real_count)+float(motif.rand_mean)!=0:
                motif_stats.add_motif(motif_id)
                motif_stats.motifs[motif_id].real = int(motif.real_count)
                motif_stats.motifs[motif_id].random_m = float(motif.rand_mean)
                motif_stats.motifs[motif_id].random_sd = float(motif.rand_std_dev)
                motif_stats.motifs[motif_id].real_z = float(motif.real_zscore)
                motif_stats.motifs[motif_id].mean_weight = 0.0
                motif_stats.motifs[motif_id].sd_weight = 0.0
                motif_stats.motifs[motif_id].median_weight = 0.0
                motif_stats.motifs[motif_id].firstq_weight = 0.0
                motif_stats.motifs[motif_id].thirdq_weight = 0.0

            motif_result = motif_result.next

    cmfinder.list64_free_mem(results)

    return motif_stats

def weighted_motif_stats(mfinderi, motif_stats, fweight, allmotifs):

    results = cmfinder.motif_participation(mfinderi)
    CI={}

    r_l = results.l
    members = cmfinder.intArray(mfinderi.MotifSize)
    while (r_l != None):
        motif = cmfinder.get_motif(r_l.p)
        id = int(motif.id)
        idx=0
        
        CI[id]=np.zeros(motif_stats.motifs[id].real)

        am_l = motif.all_members.l

        while (am_l != None):

            cmfinder.get_motif_members(am_l.p, members, mfinderi.MotifSize)
            py_members = [int(members[i]) for i in xrange(mfinderi.MotifSize)]

            weight = fweight([motif_stats.links[x].weight for x in permutations(py_members, 2) if x in motif_stats.links])

            CI[id][idx] = weight
            idx+=1

            am_l = am_l.next

        r_l = r_l.next

    cmfinder.res_tbl_mem_free_single(results)

    for motif_id in motif_stats.motifs:
        if motif_stats.motifs[motif_id].real>0:
            motif_stats.motifs[motif_id].mean_weight, motif_stats.motifs[motif_id].sd_weight, motif_stats.motifs[motif_id].median_weight, motif_stats.motifs[motif_id].thirdq_weight, motif_stats.motifs[motif_id].firstq_weight = confidence_interval(CI[motif_id])

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

    if motifsize < 2:
        sys.stderr.write("Error: this is not a valid motif size.\n")
        sys.exit()

    if motifsize > 8:
        sys.stderr.write("Error: this is not a recommended motif size.\n")
        sys.exit()

    if motifsize > 4 and allmotifs:
        sys.stderr.write("Warning: 'allmotifs' will be ignored for this motif size and motif_participation will only register existing motifs in the real network.\n")
        allmotifs=False

    # do we want to randomize the network first?
    if randomize:
        #This will restart the whole object
        network = random_network(network, usemetropolis = usemetropolis)

    # initialize the heinous input struct
    web = cmfinder.mfinder_input()

    # setup the network info
    network, stats, web.Edges, web.NumEdges = mfinder_network_setup(network)

    # add or check some basics of stats object
    if stoufferIDs and motifsize!=3:
        sys.stderr.write("Warning: 'stoufferIDs' can only be true when 'motifsize=3' in unipartite networks.\n")
        stats.stoufferIDs = False
    else:
        stats.stoufferIDs = stoufferIDs

    if stats.motifsize:
        if stats.motifsize != motifsize:
            sys.stderr.write("Error: you're trying to mix motif sizes.\n")
            sys.exit()
    else:
        stats.motifsize = motifsize

    if stats.networktype:
        stats.networktype = "unipartite"

    if stats.weighted!=None:
        if stats.weighted != weighted:
            sys.stderr.write("Warning: you're trying to mix two different motif analyses (weighted and not weighted). Be careful!\n")

    stats.weighted = weighted

    if fweight==None:
        fweight = default_fweight

    # parameterize the analysis
    web.MotifSize = motifsize
    web.Randomize = 0
    web.UseMetropolis = 0
    if len(stats.motifs) == 0:
        web.NRandomizations = 0
        web.UseMetropolis = 0
        motif_stats(web,stats, allmotifs)

    web.MaxMembersListSz = max([stats.motifs[x].real for x in stats.motifs])+1

    #TODO I can also run this inside participation
    if stats.weighted:
        weighted_motif_stats(web,stats,fweight,allmotifs)

    #check if this function has already been run
    if len(stats.nodes[stats.nodes.keys()[0]].motifs) != 0:
        for x in stats.nodes.keys():
            stats.nodes[x].motifs = dict()
        if len(stats.links[stats.links.keys()[0]].motifs) != 0:
            for x in stats.links.keys():
                stats.links[x].motifs = dict()

    return participation_stats(web,stats,links,allmotifs,fweight)


def participation_stats(mfinderi, participation, links, allmotifs, fweight):

    #TODO This is sooooo coooool!
    #possible_roles = generate_role_files(mfinderi.MotifSize, networktype=networktype, nlayers=nlayers, layers_method=layers_method)
    #possible_motifs = set([tuple([i[0], i[len(i)-1]]) for i in possible_roles])

    results = cmfinder.motif_participation(mfinderi)

    r_l = results.l
    members = cmfinder.intArray(mfinderi.MotifSize)
    while (r_l != None):
        motif = cmfinder.get_motif(r_l.p)
        id = int(motif.id)

        am_l = motif.all_members.l
        while (am_l != None):
            cmfinder.get_motif_members(am_l.p, members, mfinderi.MotifSize)
            py_members = [int(members[i]) for i in xrange(mfinderi.MotifSize)]

            if participation.weighted:

                weight = fweight([participation.links[x].weight for x in permutations(py_members, 2) if x in participation.links])

                for m in py_members:
                    try:
                        participation.nodes[m].motifs[id] += 1
                    except KeyError:
                        participation.nodes[m].motifs[id] = 1
                    try:
                        participation.nodes[m].weighted_motifs[id] += weight
                    except KeyError:
                        participation.nodes[m].weighted_motifs[id] = weight
            else:
                for m in py_members:
                    try:
                        participation.nodes[m].motifs[id] += 1
                    except KeyError:
                        participation.nodes[m].motifs[id] = 1

            if links:
                for m in permutations(py_members, 2):
                    if m in participation.links:
                        try:
                            participation.links[m].motifs[id] += 1
                        except KeyError:
                            participation.links[m].motifs[id] = 1

                        if participation.weighted:
                            try:
                                participation.links[m].weighted_motifs[id] += weight
                            except KeyError:
                                participation.links[m].weighted_motifs[id] = weight

            am_l = am_l.next

        r_l = r_l.next

    cmfinder.res_tbl_mem_free_single(results)

    possible_motifs = set(participation.motifs.keys())


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
                allroles = False,
                weighted = False,
                fweight = None,
                layers = None,
                layers_method = "complete"
                ):

    if motifsize < 2:
        sys.stderr.write("Error: this is not a valid motif size.\n")
        sys.exit()

    if (networktype=="unipartite" and motifsize>4) or (networktype=="bipartite" and motifsize>6):
        sys.stderr.write("Error: the analysis of the motif-role profiles can only be done for motif size 2, 3 and 4 in unipartite networks and up to motif size 6 in bipartite networks.\n")
        sys.exit()

    # do we want to randomize the network first?
    if randomize:
        #This will restart the whole object
        network = random_network(network, usemetropolis = usemetropolis)

    # initialize the heinous input struct
    web = cmfinder.mfinder_input()

    # setup the network info
    network, stats, web.Edges, web.NumEdges = mfinder_network_setup(network)

    # add or check some basics of stats object
    if stoufferIDs and motifsize!=3:
        sys.stderr.write("Warning: 'stoufferIDs' can only be true when 'motifsize=3' in unipartite networks.\n")
        stats.stoufferIDs = False
    else:
        stats.stoufferIDs = stoufferIDs

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

    if stats.weighted!=None:
        if stats.weighted != weighted:
            sys.stderr.write("You're trying to mix two different motif analyses (weighted and not weighted). Be careful!\n")

    stats.weighted = weighted

    if fweight==None:
        fweight = default_fweight

    if layers!=None:
        if layers_method not in ["complete", "simple"]:
            sys.stderr.write("Error: the variable layers_method should be set to either 'complete' or 'simple'. Default 'complete'.\n")
            sys.exit()
        nlayers, stats = read_layers(layers, stats)
    else:
        nlayers=1

    # parameterize the analysis
    web.MotifSize = motifsize
    web.Randomize = 0
    web.UseMetropolis = 0
    if len(stats.motifs) == 0:
        web.NRandomizations = 0
        web.UseMetropolis = 0
        motif_stats(web,stats, allroles)

    web.MaxMembersListSz = max([stats.motifs[x].real for x in stats.motifs])+1

    #TODO I can also run this inside role_stats
    if stats.weighted:
        weighted_motif_stats(web,stats,fweight,allroles)

    #check if this function has already been run. TODO this breaks...
    if len(stats.nodes[stats.nodes.keys()[0]].roles) != 0:
        for x in stats.nodes.keys():
            stats.nodes[x].roles = dict()
        if len(stats.links[stats.links.keys()[0]].roles) != 0:
            for x in stats.links.keys():
                stats.links[x].roles = dict()

    # determine all nodes' role statistics
    return role_stats(web,stats,links,networktype,allroles,fweight,nlayers,layers_method)


def role_stats(mfinderi,roles,links,networktype,allroles,fweight, nlayers,layers_method):
    results = cmfinder.motif_participation(mfinderi)
    actual_roles = set([])

    possible_roles = generate_role_files(mfinderi.MotifSize, networktype=networktype, nlayers=nlayers, layers_method=layers_method)


    if links:
        possible_linkroles = set([])
        actual_linkroles = set([])
        if networktype == "unipartite":
            for m,l in UNIPARTITE_LINKS_ROLES[mfinderi.MotifSize]:
                possible_linkroles.update([tuple([m] + list(x)) for x in l])
        elif networktype == "bipartite":
            for m,l in BIPARTITE_LINKS_ROLES[mfinderi.MotifSize]:
                possible_linkroles.update([tuple([m] + list(x)) for x in l])

    r_l = results.l
    members = cmfinder.intArray(mfinderi.MotifSize)
    while (r_l != None):
        motif = cmfinder.get_motif(r_l.p)
        id = int(motif.id)
        
        am_l = motif.all_members.l
        while (am_l != None):
            cmfinder.get_motif_members(am_l.p, members, mfinderi.MotifSize)
            py_members = [int(members[i]) for i in xrange(mfinderi.MotifSize)]

            py_motif = [x for x in permutations(py_members, 2) if x in roles.links]

            if roles.weighted:
                weight_motif = fweight([roles.links[x].weight for x in py_motif])
                weight_i = [fweight([roles.links[x].weight for x in py_motif if x[0]==m or x[1]==m]) for m in py_members]
                weight = weight_motif/float(sum(weight_i))

            if nlayers==1:
                for idm, m in enumerate(py_members):
                    npred, nprey = 0, 0
                    for othernode in py_members:
                        if (othernode,m) in py_motif:
                            npred+=1
                        if (m,othernode) in py_motif:
                            nprey+=1
                    key = (id, npred, nprey)

                    if key not in possible_roles:
                        connected_to = set([othernode for othernode in py_members if othernode != m and (othernode,m) in py_motif])
                        npreys = np.sort([sum([(i,j) in py_motif for j in py_members if j != i]) for i in connected_to])
                        connected_to = set([othernode for othernode in py_members if othernode != m and (m,othernode) in py_motif])
                        npreds = np.sort([sum([(j,i) in py_motif for j in py_members if j != i]) for i in connected_to])
                        key = tuple([id, tuple(npreys), tuple(npreds)])

                    if key not in possible_roles:
                        print >> sys.stderr, key
                        print >> sys.stderr, "Apparently there is a role you aren't accounting for in 'roles.py'. "
                        sys.exit()

                    try:
                        roles.nodes[m].roles[key] += 1
                    except KeyError:
                        roles.nodes[m].roles[key] = 1
    
                    if roles.weighted:
                        try:
                            roles.nodes[m].weighted_roles[key] += weight_i[idm]*weight
                        except KeyError:
                            roles.nodes[m].weighted_roles[key] = weight_i[idm]*weight

                    actual_roles.add(key)

            else:
                key=dict()
                extra=dict()
                for idm, m in enumerate(py_members):
                    key[idm] = [id, roles.nodes[m].layer]
                    for ly in range(1,nlayers+1):
                        npred, nprey = 0, 0
                        for othernode in py_members:
                            if (othernode,m) in py_motif and roles.nodes[othernode].layer==ly:
                                npred+=1
                            if (m,othernode) in py_motif and roles.nodes[othernode].layer==ly:
                                nprey+=1
                        key[idm] += [npred, nprey]
                    key[idm] = tuple(key[idm])


                    extra[idm] = [id, roles.nodes[m].layer]
                    for ly in range(1,nlayers+1):
                        connected_to = set([othernode for othernode in py_members if ((othernode != m) and ((othernode,m) in py_motif) and (roles.nodes[othernode].layer==ly))])
                        npreys = np.sort([sum([(i,j) in py_motif for j in py_members if j != i]) for i in connected_to])
                        connected_to = set([othernode for othernode in py_members if ((othernode != m) and ((m,othernode) in py_motif) and (roles.nodes[othernode].layer==ly))])
                        npreds = np.sort([sum([(j,i) in py_motif for j in py_members if j != i]) for i in connected_to])
                        extra[idm] += [tuple(npreys), tuple(npreds)]
                    extra[idm] = tuple(extra[idm])
                
                lkey = int("".join([str(_i[0]) for _i in sorted(extra.values())]))

                for idm, m in enumerate(py_members):
                    if key[m] not in possible_roles:
                        _key = tuple(list(key[m])+[lkey])
                    else:
                        _key = tuple(list(extra[m])+[lkey])

                    if _key not in possible_roles:
                        print >> sys.stderr, _key
                        print >> sys.stderr, "Apparently there is a role you aren't accounting for in 'roles.py'. "
                        sys.exit()

                    try:
                        roles.nodes[m].roles[_key] += 1
                    except KeyError:
                        roles.nodes[m].roles[_key] = 1
    
                    if roles.weighted:
                        try:
                            roles.nodes[m].weighted_roles[_key] += weight_i[idm]*weight
                        except KeyError:
                            roles.nodes[m].weighted_roles[_key] = weight_i[idm]*weight

                    actual_roles.add(_key)

            if links:
                for idm, m in enumerate(py_members):
                    for n in py_members:
                        if n!=m and (n,m) in roles.links: 
                            npred1,npred2,nprey1,nprey2 = 0,0,0,0
                            for othernode in py_members:
                                if (othernode,n) in py_motif:
                                    npred1+=1
                                if (n,othernode) in py_motif:
                                    nprey1+=1
                                if (othernode,m) in py_motif:
                                    npred2+=1
                                if (m,othernode) in py_motif:
                                    nprey2+=1

                            key = (id, (npred1, nprey1),(npred2, nprey2))

                            if key not in possible_linkroles:
                                key = (id, (npred2, nprey2),(npred1, nprey1))

                            # There are three motifs containing links that cannot be uniquely specified by (npred1,nprey1),(npred2,nprey2).
                            if key not in possible_linkroles:
                                nconnected=set([othernode for othernode in py_members if othernode != n and (n,othernode) in py_motif])
                                mconnected=set([othernode for othernode in py_members if othernode != m and (othernode,m) in py_motif])

                                npreypreds=sorted([sum([(i,j) in py_motif for i in py_members if i!=j]) for j in nconnected]) # predators for each prey of n
                                mpredpreys=sorted([sum([(i,j) in py_motif for j in py_members if j!=i]) for i in mconnected]) # prey for each predator of m

                                # One link has both pred and prey with nonunique roles. 
                                key = (id, (npred1, nprey1,tuple(npreypreds)),(npred2,nprey2,tuple(mpredpreys)))
                                # One link has only non-unique predator
                                if key not in possible_linkroles:
                                    key= (id, (npred1, nprey1,tuple(npreypreds)),(npred2,nprey2))
                                # One link has only non-unique prey
                                if key not in possible_linkroles:
                                    key = (id, (npred1, nprey1),(npred2,nprey2,tuple(mpredpreys)))

                            if key not in possible_linkroles:
                                print >> sys.stderr, key
                                print >> sys.stderr, "Apparently there is a link you aren't accounting for in roles.py."

                            try:
                                roles.links[(n,m)].roles[key] += 1
                            except KeyError:
                                roles.links[(n,m)].roles[key] = 1

                            if roles.weighted:
                                try:
                                    roles.links[(n,m)].weighted_roles[key] += weight_motif
                                except KeyError:
                                    roles.links[(n,m)].weighted_roles[key] = weight_motif

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

            if roles.weighted:
                try:
                    x = roles.nodes[n].weighted_roles[r]
                except KeyError:
                    roles.nodes[n].weighted_roles[r] = 0

    if links:
        for n in roles.links:
            for r in possible_linkroles:
                try:
                    x = roles.links[n].roles[r]
                except KeyError:
                    roles.links[n].roles[r] = 0

                if roles.weighted:
                    try:
                        x = roles.links[n].weighted_roles[r]
                    except KeyError:
                        roles.links[n].weighted_roles[r] = 0

    return roles



##############################################################
##############################################################
# RUN IT ALL!
##############################################################
##############################################################

def pymfinder(network,
              links=False,
              motifsize = 3,
              stoufferIDs = False,
              allmotifs = False,
              nrandomizations = 0,
              randomize = False,
              usemetropolis = False,
              networktype = "unipartite",
              weighted = False
              ):


    stats = motif_participation(network,
                                links = links,
                                motifsize = motifsize,
                                randomize = randomize,
                                usemetropolis = usemetropolis,
                                stoufferIDs = stoufferIDs,
                                allmotifs = allmotifs,
                                weighted = weighted)

    stats = motif_roles(stats,
                        links = links,
                        motifsize = motifsize,
                        randomize = False,
                        usemetropolis = usemetropolis,
                        stoufferIDs = stoufferIDs,
                        networktype = networktype,
                        allroles = allmotifs,
                        weighted = weighted)

    if nrandomizations != 0:
        stats = motif_structure(stats,
                                motifsize = motifsize,
                                nrandomizations = nrandomizations,
                                usemetropolis = usemetropolis,
                                stoufferIDs = stoufferIDs,
                                weighted = weighted)



    return stats

##############################################################
##############################################################
# C'est fini
##############################################################
##############################################################
