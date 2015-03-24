
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
