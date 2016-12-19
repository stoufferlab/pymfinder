
import unittest

from pymfinder import motif_structure
from pymfinder import motif_participation
from pymfinder import motif_roles

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
		motifs=[result[n]["real"] for n in result]
		self.assertTrue(motifs==[1]*len(motifs))

    def test_unipartite_motif_participation(self):
	for i in range(0,2):
		result = motif_participation(self.test_filename_u[i], motifsize=i+2,)
		all_roles = result[result.keys()[0]].keys()
		roles=[[result[n][r] for r in all_roles] for n in result]
		check_columns=[sum(x) for x in zip(*roles)]==[i+2]*len(all_roles)
		check_rows=[sum(x) for x in roles]==[1]*len(result)
        	self.assertTrue(check_columns and check_rows)

    def test_unipartite_motif_roles(self):
	for i in range(0,2):
		result = motif_roles(self.test_filename_u[i], motifsize=i+2,)
		all_roles = result[result.keys()[0]].keys()
		roles=[[result[n][r] for r in all_roles] for n in result]
		columns=[sum(x) for x in zip(*roles)]
		check_total_nodes=sum(columns)==3*13 or 2*2
		check_columns=any(v==0 for v in columns)==False
		check_rows=[sum(x) for x in roles]==[1]*len(result)
        	self.assertTrue(check_columns and check_rows and check_total_nodes)

    def test_bipartite_motif_roles(self):
	for i in range(0,5):
		result = motif_roles(self.test_filename_b[i], motifsize=i+2, networktype = "bipartite",)
		all_roles = result[result.keys()[0]].keys()
		roles=[[result[n][r] for r in all_roles] for n in result]
		check_columns=any(v==0 for v in [ sum(x) for x in zip(*roles)])==False
		check_rows=[sum(x) for x in roles]==[1]*len(result)
        	self.assertTrue(check_columns and check_rows)

def suite():
    suite = unittest.TestLoader().loadTestsFromTestCase(pymfinderTestCase)
    return suite

if __name__ == '__main__':
    unittest.TextTestRunner(verbosity=2).run(suite())
