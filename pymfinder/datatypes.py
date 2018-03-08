
import sys

##############################################################
##############################################################
# What is a motif?
##############################################################
##############################################################

class Motif(object):
    """Motif info class"""

    def __init__(self,motif_id=None):
        self.id = motif_id
        self.real = None
        self.all_members = set()
        self.random = list()
        self.random_m = None
        self.random_sd = None
        self.real_z = None
        self.weighted = None

    # def __rep__(self):
    #     return(str(self.id))


##############################################################
##############################################################
# What is a node or a link?
##############################################################
##############################################################

class NodeLink(object):
    """NodeLink info class"""

    def __init__(self,nodelink_id=None):
        self.id = nodelink_id
        self.motifs = dict()
        self.roles = dict()
        self.weighted_roles = dict()
        self.weighted_motifs = dict()


##############################################################
##############################################################
# What is a network stats?
##############################################################
##############################################################

class NetworkStats(object):
    """NetworkStats summary info class"""

    def __init__(self, motifsize = None, networktype = None, weighted = None):
        self.motifs = dict()
        self.nodes = dict()
        self.links = dict()
        self.networktype = networktype
        self.motifsize = motifsize
        self.weighted = weighted

    def add_motif(self,motif_id):
        if motif_id in self.motifs:
            sys.stderr.write("You're trying to add a motif more than once. According to the developers, this is classified as highly unusual.\n")
        else:
            self.motifs[motif_id] = Motif(motif_id)

    def add_node(self, node_id, node_name=None):
        if node_id in self.nodes:
            sys.stderr.write("You're trying to add a node more than once. According to the developers, this is classified as highly unusual.\n")
        else:
            self.nodes[node_id] = NodeLink(node_name)

    def add_link(self, link_id, link_name=None):
        if link_id in self.links:
            sys.stderr.write("You're trying to add a link more than once. According to the developers, this is classified as highly unusual.\n")
        else:
            self.links[link_id] = NodeLink(link_name)

    def use_stouffer_IDs(self):
        from roles import STOUFFER_MOTIF_IDS, UNIPARTITE_ROLES, UNIPARTITE_LINKS_ROLES

        ineligible_ids = [motif_id for motif_id in self.motifs if motif_id not in STOUFFER_MOTIF_IDS]

        if len(ineligible_ids) == 0:
            self.motifs = dict([(STOUFFER_MOTIF_IDS[motif_id],self.motifs[motif_id]) for motif_id in self.motifs])
        else:
            pass

        if self.motifsize == 3 and self.networktype == "unipartite":
            possible_roles = []
            for motif,roles in UNIPARTITE_ROLES[self.motifsize]:
                possible_roles += [tuple([motif] + list(role)) for role in roles]

            for n in self.nodes:
                self.nodes[n].motifs = dict([(STOUFFER_MOTIF_IDS[id],self.nodes[n].motifs[id]) for id in self.nodes[n].motifs])
                self.nodes[n].roles = dict([(possible_roles.index(id)+1,self.nodes[n].roles[id]) for id in self.nodes[n].roles])
        else:
            pass

        if self.motifsize == 3 and self.networktype == "unipartite":
            possible_roles = []
            for m,r in UNIPARTITE_LINKS_ROLES[self.motifsize]:
                possible_roles += [tuple([m] + list(x)) for x in r]

            for n in self.links:
                self.links[n].motifs = dict([(STOUFFER_MOTIF_IDS[id],self.links[n].motifs[id]) for id in self.links[n].motifs])
                self.links[n].roles = dict([(possible_roles.index(id)+1,self.links[n].roles[id]) for id in self.links[n].roles])
        else:
            pass

    # DEBUG: it would be nice to be able to turn the header on and off
    def __str__(self):

        output = ""

        if self.motifs != dict():
            output = output + " ".join(['motif',
                               'real',
                               'rand',
                               'srand',
                               'zscore',
                               'weighted',]) + '\n'

            # set up the data itself
            for m in sorted(self.motifs.keys()):
                output = output + " ".join(["%s" % str(m),
                                   "%i" % self.motifs[m].real,
                                   "%.3f" % self.motifs[m].random_m,
                                   "%.3f" % self.motifs[m].random_sd,
                                   "%.3f" % self.motifs[m].real_z,
                                   "%.3f" % self.motifs[m].weighted,
                                  ]) + '\n'
            output = output + '\n'

        if self.nodes[self.nodes.keys()[0]].motifs != dict():
            # set up a header
            output = output + " ".join(["node"]+list(map(str,sorted(self.nodes[self.nodes.keys()[0]].motifs.keys()))))
            output = output + '\n'

            # set up the data itself
            # set up a header
            if self.weighted:
                for m in sorted(self.nodes.keys()):
                    output = output + " ".join([str(self.nodes[m].id)] + list(map(str,[j for i,j in sorted(self.nodes[m].weighted_motifs.items())]))) + '\n'
                output = output + '\n'
            else:
                for m in sorted(self.nodes.keys()):
                    output = output + " ".join([str(self.nodes[m].id)] + list(map(str,[j for i,j in sorted(self.nodes[m].motifs.items())]))) + '\n'
                output = output + '\n'

            if self.links[self.links.keys()[0]].motifs != dict():
                # set up a header
                output = output + " ".join(["link"]+list(map(str,sorted(self.links[self.links.keys()[0]].motifs.keys()))))
                output = output + '\n'

                if self.weighted:
                    # set up the data itself
                    for m in sorted(self.links.keys()):
                        output = output + " ".join([str(self.links[m].id)] + list(map(str,[j for i,j in sorted(self.links[m].weighted_motifs.items())]))) + '\n'
                    output = output + '\n'
                else:
                    # set up the data itself
                    for m in sorted(self.links.keys()):
                        output = output + " ".join([str(self.links[m].id)] + list(map(str,[j for i,j in sorted(self.links[m].motifs.items())]))) + '\n'
                    output = output + '\n'


        if self.nodes[self.nodes.keys()[0]].roles != dict():
            # set up a header
            # DEBUG: consider changing role for ".".join(map(str, role)) fixing problems when STOUFFERID=True
            output = output+" ".join(["node"]+list(map(str,[role for role in sorted(self.nodes[self.nodes.keys()[0]].roles.keys())])))
            output = output + '\n'

            # set up the data itself
            for m in sorted(self.nodes.keys()):
                output = output + " ".join([str(m)] + list(map(str,[j for i,j in sorted(self.nodes[m].roles.items())]))) + '\n'
            output = output + '\n'

            if self.links[self.links.keys()[0]].roles != dict():
                # set up a header
                output = output + " ".join(["link"]+list(map(str,sorted(self.links[self.links.keys()[0]].roles.keys()))))
                output = output + '\n'

                # set up the data itself
                for m in sorted(self.links.keys()):
                    output = output + " ".join([str(m)] + list(map(str,[j for i,j in sorted(self.links[m].roles.items())]))) + '\n'
                output = output + '\n'

        #TODO Clean all the print functions... it's ugly AF

        # return this ghastly beast
        return(output)


