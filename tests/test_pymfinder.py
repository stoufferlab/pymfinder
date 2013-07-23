
import unittest

from pymfinder import motif_structure
from pymfinder import motif_participation
from pymfinder import motif_roles

import test_results

class pymfinderTestCase(unittest.TestCase):
    def setUp(self):
        import os
        self.test_filename = os.path.dirname(__file__) + "/unipartite-test.net"
        try:
            with open(self.test_filename) as file:
                pass
        except IOError:
            import sys
            sys.stderr.write("Cannot find the test network 'unipartite-test.net'.\n")
            sys.exit(1)

    def test_motif_structure(self):
        result = motif_structure(self.test_filename,motifsize=3,nrandomizations=0)
        self.assertEqual(result, test_results.motif_structure)
    
    def test_motif_participation(self):
        result = motif_participation(self.test_filename,motifsize=3,stoufferIDs=True)
        self.assertEqual(result, test_results.motif_participation)

    def test_motif_roles(self):
        result = motif_roles(self.test_filename,motifsize=3,stoufferIDs=True)
        self.assertEqual(result, test_results.motif_roles)

def suite():
    suite = unittest.TestLoader().loadTestsFromTestCase(pymfinderTestCase)
    return suite

if __name__ == '__main__':
    unittest.TextTestRunner(verbosity=2).run(suite())
