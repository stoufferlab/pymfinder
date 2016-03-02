
import unittest

from pymfinder import motif_structure
from pymfinder import motif_participation
from pymfinder import motif_roles

import test_results

class pymfinderTestCase(unittest.TestCase):
    def setUp(self):
        import os
        self.test_filename1 = os.path.dirname(__file__) + "/../pymfinder/data/unipartite-test.net"
        self.test_filename2 = [os.path.dirname(__file__) + "/../pymfinder/data/bipartite-"+str(x)+"-test.net" for x in range(2,7)]

        try:
            with open(self.test_filename1) as file:
                pass
        except IOError:
            import sys
            sys.stderr.write("Cannot find the test network 'unipartite-test.net'.\n")
            sys.exit(1)

	for x in range(0,5):
		try:
		    with open(self.test_filename2[x]) as file:
		        pass
		except IOError:
		    import sys
		    sys.stderr.write("Cannot find the test network 'bipartite-"+str(x)+"-test.net'.\n")
		    sys.exit(1)

    def test_motif_structure(self):
        result = motif_structure(self.test_filename1,motifsize=3,nrandomizations=0)
        self.assertEqual(result, test_results.motif_structure)
    
    def test_motif_participation(self):
        result = motif_participation(self.test_filename1,motifsize=3,stoufferIDs=True)
        self.assertEqual(result, test_results.motif_participation)

    def test_motif_roles(self):
        result = motif_roles(self.test_filename1,motifsize=3,stoufferIDs=True)
        self.assertEqual(result, test_results.motif_roles)

    def test_bipartite_motif_roles(self):
	for i in range(0,5):
		result = motif_roles(self.test_filename2[i], motifsize=i+2, networktype = "bipartite",)
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
