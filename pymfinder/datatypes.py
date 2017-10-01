
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

    # def __rep__(self):
    #     return(str(self.id))

##############################################################
##############################################################
# What is a motif stat?
##############################################################
##############################################################

class MotifStats(object):
    """MotifStats summary info class"""

    def __init__(self):
        self.motifs = dict()

    def add_motif(self,motif_id):
        if motif_id in self.motifs:
            sys.stderr.write("You're trying to add a motif more than once. According to the developers, this is classified as highly unusual.\n")
        else:
            self.motifs[motif_id] = Motif(motif_id)

    def use_stouffer_IDs(self):
        from roles import STOUFFER_MOTIF_IDS
        ineligible_ids = [motif_id for motif_id in self.motifs if motif_id not in STOUFFER_MOTIF_IDS]

        if len(ineligible_ids) == 0:
            self.motifs = dict([(STOUFFER_MOTIF_IDS[motif_id],self.motifs[motif_id]) for motif_id in self.motifs])
        else:
            pass

    # DEBUG: it would be nice to be able to turn the header on and off
    def __str__(self):
        # set up a header
        output = " ".join(['motif',
                           'real',
                           'rand',
                           'srand',
                           'zscore',]) + '\n'

        # set up the data itself
        for m in sorted(self.motifs.keys()):
            output = output + " ".join(["%s" % str(m),
                               "%i" % self.motifs[m].real,
                               "%.3f" % self.motifs[m].random_m,
                               "%.3f" % self.motifs[m].random_sd,
                               "%.3f" % self.motifs[m].real_z,
                              ]) + '\n'

        # return this ghastly beast
        return(output)


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


##############################################################
##############################################################
# What is a node stat?
##############################################################
##############################################################

class NodeStats(object):
    """ParticipationStats summary info class"""

    def __init__(self, motifsize = 3, networktype = "unipartite"):
        self.nodes = dict()
	self.links = dict()
        self.networktype = networktype
        self.motifsize = motifsize

    def add_node(self,node_id):
        if node_id in self.nodes:
            sys.stderr.write("You're trying to add a node more than once. According to the developers, this is classified as highly unusual.\n")
        else:
            self.nodes[node_id] = NodeLink(node_id)

    def add_link(self,link_id):
        if link_id in self.links:
            sys.stderr.write("You're trying to add a link more than once. According to the developers, this is classified as highly unusual.\n")
        else:
            self.links[link_id] = NodeLink(link_id)

    def use_stouffer_IDs(self):
        from roles import STOUFFER_MOTIF_IDS, UNIPARTITE_ROLES

        possible_roles = []
        for motif,roles in UNIPARTITE_ROLES[self.motifsize]:
            possible_roles += [tuple([motif] + list(role)) for role in roles]

        for n in self.nodes:
            if self.motifsize == 3 and self.networktype == "unipartite":
                self.nodes[n].motifs = dict([(STOUFFER_MOTIF_IDS[id],self.nodes[n].motifs[id]) for id in self.nodes[n].motifs])
                self.nodes[n].roles = dict([(possible_roles.index(id)+1,self.nodes[n].roles[id]) for id in self.nodes[n].roles])
            else:
                pass

        for n in self.links:
            if self.motifsize == 3 and self.networktype == "unipartite":
                self.links[n].motifs = dict([(STOUFFER_MOTIF_IDS[id],self.links[n].motifs[id]) for id in self.links[n].motifs])
                self.links[n].roles = dict([(possible_roles.index(id)+1,self.links[n].roles[id]) for id in self.links[n].roles])
            else:
                pass

    # DEBUG: it would be nice to be able to turn the header on and off
    def __str__(self):

        if self.nodes[self.nodes.keys()[0]].motifs != dict():
            # set up a header
            output = " ".join(["node"]+list(map(str,sorted(self.nodes[self.nodes.keys()[0]].motifs.keys()))))
            output = output + '\n'

            # set up the data itself
            for m in sorted(self.nodes.keys()):
                output = output + " ".join([str(m)] + list(map(str,[j for i,j in sorted(self.nodes[m].motifs.items())]))) + '\n'

            if self.links[self.links.keys()[0]].motifs != dict():
                # set up a header
                output = output + '\n\n' + " ".join(["link"]+list(map(str,sorted(self.links[self.links.keys()[0]].motifs.keys()))))
                output = output + '\n'

                # set up the data itself
                for m in sorted(self.links.keys()):
                    output = output + " ".join([str(m)] + list(map(str,[j for i,j in sorted(self.links[m].motifs.items())]))) + '\n'

        if self.nodes[self.nodes.keys()[0]].roles != dict():
            # set up a header
            # DEBUG: consider changing role for ".".join(map(str, role)) fixing problems when STOUFFERID=True
            output = " ".join(["node"]+list(map(str,[role for role in sorted(self.nodes[self.nodes.keys()[0]].roles.keys())])))
            output = output + '\n'

            # set up the data itself
            for m in sorted(self.nodes.keys()):
                output = output + " ".join([str(m)] + list(map(str,[j for i,j in sorted(self.nodes[m].roles.items())]))) + '\n'

            if self.links[self.links.keys()[0]].roles != dict():
                # set up a header
                output = output + '\n\n' + " ".join(["link"]+list(map(str,sorted(self.links[self.links.keys()[0]].roles.keys()))))
                output = output + '\n'

                # set up the data itself
                for m in sorted(self.links.keys()):
                    output = output + " ".join([str(m)] + list(map(str,[j for i,j in sorted(self.links[m].roles.items())]))) + '\n'

        # return this ghastly beast
        return(output)


