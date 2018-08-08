
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
        self.mean_weight = None
        self.sd_weight = None
        self.median_weight = None
        self.thirdq_weight = None
        self.firstq_weight = None

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

    def __init__(self, motifsize = None, networktype = None, weighted = None, stoufferIDs = None):
        self.motifs = dict()
        self.nodes = dict()
        self.links = dict()
        self.networktype = networktype
        self.motifsize = motifsize
        self.weighted = weighted
        self.stoufferIDs = stoufferIDs

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

    # DEBUG: it would be nice to be able to turn the header on and off
    def __str__(self):

        from roles import STOUFFER_MOTIF_IDS

        if self.stoufferIDs:
            ineligible_ids = [motif_id for motif_id in self.motifs if motif_id not in STOUFFER_MOTIF_IDS]
        else:
            ineligible_ids = [1,2,3]

        output = ""

        if self.motifs != dict():
            output = output + " ".join(['motif',
                               'real',
                               'rand',
                               'srand',
                               'zscore',
                               'weight-mean',
                               'weight-sd',]) + '\n'

            # set up the data itself
            for m in sorted(self.motifs.keys()):

                if len(ineligible_ids) == 0:
                    motifname = STOUFFER_MOTIF_IDS[m]
                else:
                    motifname = m

                output = output + " ".join(["%s" % str(motifname),
                                   "%i" % self.motifs[m].real,
                                   "%.3f" % self.motifs[m].random_m,
                                   "%.3f" % self.motifs[m].random_sd,
                                   "%.3f" % self.motifs[m].real_z,
                                   "%.3f" % self.motifs[m].mean_weight,
                                   "%.3f" % self.motifs[m].sd_weight,
                                  ]) + '\n'
            output = output + '\n'

        if self.nodes[self.nodes.keys()[0]].motifs != dict():
            # set up a header
            if len(ineligible_ids) == 0:
                output = output + " ".join(["node"]+list(map(str, [STOUFFER_MOTIF_IDS[mid] for mid in sorted(self.nodes[self.nodes.keys()[0]].motifs.keys())])))
            else:
                output = output + " ".join(["node"]+list(map(str,sorted(self.nodes[self.nodes.keys()[0]].motifs.keys()))))

            output = output + '\n'

            # set up the data itself
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
            output = output+" ".join(["node"]+list(map(str,[role for role in sorted(self.nodes[self.nodes.keys()[0]].roles.keys())])))
            output = output + '\n'

            if self.weighted:
                for m in sorted(self.nodes.keys()):
                    output = output + " ".join([str(self.nodes[m].id)] + list(map(str,[j for i,j in sorted(self.nodes[m].weighted_roles.items())]))) + '\n'
                output = output + '\n'
            else:
                for m in sorted(self.nodes.keys()):
                    output = output + " ".join([str(self.nodes[m].id)] + list(map(str,[j for i,j in sorted(self.nodes[m].roles.items())]))) + '\n'
                output = output + '\n'


            if self.links[self.links.keys()[0]].roles != dict():
                # set up a header
                output = output + " ".join(["link"]+list(map(str,sorted(self.links[self.links.keys()[0]].roles.keys()))))
                output = output + '\n'

                if self.weighted:
                    # set up the data itself
                    for m in sorted(self.links.keys()):
                        output = output + " ".join([str(self.links[m].id)] + list(map(str,[j for i,j in sorted(self.links[m].weighted_roles.items())]))) + '\n'
                    output = output + '\n'
                else:
                    # set up the data itself
                    for m in sorted(self.links.keys()):
                        output = output + " ".join([str(self.links[m].id)] + list(map(str,[j for i,j in sorted(self.links[m].roles.items())]))) + '\n'
                    output = output + '\n'


        # return this ghastly beast
        return(output)


