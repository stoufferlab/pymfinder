
import unittest

from pymfinder import motif_structure
from pymfinder import motif_participation
from pymfinder import motif_roles

links_per_motif={
                 2:{2:1,6:2}, # One-way motif, two-way motif
                 3:{12:2,38:3,98:3,36:2,6:2, # One-way motifs only
                 46:4,108:4,14:3,74:3,102:4,238:6,110:5,78:4} #Motifs with 2-way ints
                 }

class pymfinderTestCase(unittest.TestCase):
    def setUp(self):
        import sys
        import os
        self.test_filename_u = [os.path.dirname(__file__) + "/../pymfinder/data/unipartite-"+str(x)+"-test.net" for x in range(2,4)]
        self.test_filename_b = [os.path.dirname(__file__) + "/../pymfinder/data/bipartite-"+str(x)+"-test.net" for x in range(2,7)]

        for x in range(0,2):
            print self.test_filename_u[x]
            try:
                with open(self.test_filename_u[x]) as file:
                    pass
            except IOError:
                sys.stderr.write("Cannot find the test network 'unipartite-"+str(x+2)+"-test.net'.\n")
                sys.exit(1)

        for x in range(0,5):
            print self.test_filename_b[x]
            try:
                with open(self.test_filename_b[x]) as file:
                    pass
            except IOError:
                sys.stderr.write("Cannot find the test network 'bipartite-"+str(x+2)+"-test.net'.\n")
                sys.exit(1)

    def test_unipartite_motif_structure(self):
        for i in range(0,2):
            result = motif_structure(self.test_filename_u[i], motifsize=i+2, nrandomizations=0)
            motifs=[result.motifs[n].real for n in result.motifs]
            self.assertTrue(motifs==[1]*len(motifs))

    def test_unipartite_motif_weighted_structure(self):
        for i in range(0,2):
            result = motif_structure(self.test_filename_u[i], motifsize=i+2, nrandomizations=0, weighted=True)
            motifs=[int(result.motifs[n].mean_weight)==n for n in result.motifs.keys()]
            self.assertTrue(all(motifs))

    def test_unipartite_motif_participation(self):
        for i in range(0,2):
            result = motif_participation(self.test_filename_u[i], motifsize=i+2,)
            all_roles = result.nodes[result.nodes.keys()[0]].motifs.keys()
            roles=[[result.nodes[n].motifs[r] for r in all_roles] for n in result.nodes]
            check_columns=[sum(x) for x in zip(*roles)]==[i+2]*len(all_roles) # Check the number of nodes in each motif
            check_rows=[sum(x) for x in roles]==[1]*len(result.nodes) # Check that each node is in only one motif
            self.assertTrue(check_columns and check_rows)

    def test_unipartite_motif_weighted_participation(self):
        for i in range(0,2):
            result = motif_participation(self.test_filename_u[i], motifsize=i+2, weighted=True)
            all_roles = result.nodes[result.nodes.keys()[0]].weighted_motifs.keys()
            roles=[[int(result.nodes[n].weighted_motifs[r]) for r in all_roles] for n in result.nodes]
            check_columns=[sum(x)/(i+2) for x in zip(*roles)]==all_roles # Check the number of nodes in each motif
            check_rows=[sum(x) for x in roles]==[int(result.nodes[x].id[:-1]) for x in result.nodes] # Check that each node is in only one motif
            self.assertTrue(check_columns and check_rows)

    def test_unipartite_link_participation(self):
        for i in range(0,2):
            result = motif_participation(self.test_filename_u[i], motifsize=i+2,links=True)
            all_roles = result.links[result.links.keys()[0]].motifs.keys() # List of motifs
            roles=[[result.links[n].motifs[r] for r in sorted(all_roles)] for n in result.links] # r are motifs, n are link tuples
            motiflist=sorted(result.motifs.keys()) # Need to keep these in the same order
            check_columns=[sum(x) for x in zip(*roles)]==[links_per_motif[i+2][motif] for motif in motiflist] # Each motif contains the number of links I think it should
            check_rows=[sum(x) for x in roles]==[1]*len(result.links) # Each link is in only one motif
            self.assertTrue(check_columns and check_rows)

    def test_unipartite_motif_roles(self):
        for i in range(0,2):
            result = motif_roles(self.test_filename_u[i], motifsize=i+2,)
            all_roles = result.nodes[result.nodes.keys()[0]].roles.keys()
            roles=[[result.nodes[n].roles[r] for r in all_roles] for n in result.nodes]
            columns=[sum(x) for x in zip(*roles)]
            check_total_nodes=sum(columns)==3*13 or 2*2
            check_columns=any(v==0 for v in columns)==False
            check_rows=[sum(x) for x in roles]==[1]*len(result.nodes)
            self.assertTrue(check_columns and check_rows and check_total_nodes)

    def test_unipartite_link_roles(self):
        for i in range(0,2):
            result = motif_roles(self.test_filename_u[i], motifsize=i+2,links=True)
            all_roles = result.links[result.links.keys()[0]].roles.keys()
            roles=[[result.links[n].roles[r] for r in all_roles] for n in result.links]
            columns=[sum(x) for x in zip(*roles)]
            #Are there the number of links we think there are?
            check_total_links=sum(columns)==sum(links_per_motif[i+2].values())
            check_columns=any(v==0 for v in columns)==False # Every link has a role?
            check_rows=[sum(x) for x in roles]==[1]*len(result.links) # Every link is in one position only
            self.assertTrue(check_columns and check_rows and check_total_links)

    def test_bipartite_motif_roles(self):
        for i in range(0,5):
            result = motif_roles(self.test_filename_b[i], motifsize=i+2, networktype = "bipartite",)
            all_roles = result.nodes[result.nodes.keys()[0]].roles.keys()
            roles=[[result.nodes[n].roles[r] for r in all_roles] for n in result.nodes]
            check_columns=any(v==0 for v in [ sum(x) for x in zip(*roles)])==False
            check_rows=[sum(x) for x in roles]==[1]*len(result.nodes)
            self.assertTrue(check_columns and check_rows)

# Bipartite link roles are not yet defined.


def suite():
    suite = unittest.TestLoader().loadTestsFromTestCase(pymfinderTestCase)
    return suite

if __name__ == '__main__':
    unittest.TextTestRunner(verbosity=2).run(suite())
