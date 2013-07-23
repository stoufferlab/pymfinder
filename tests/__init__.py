
def suite():
	import unittest
	suite = unittest.TestSuite()

	import test_pymfinder
	suite.addTests(test_pymfinder.suite())
	
	#import doctest
	#suite.addTests(doctest.DocTestSuit(motif_structure))
	#suite.addTests(doctest.DocTestSuit(motif_participation))
	#suite.addTests(doctest.DocTestSuit(motif_roles))
	
	return suite

if __name__ == '__main__':
	unittest.TextTestRunner(verbosity=2).run(suite())
