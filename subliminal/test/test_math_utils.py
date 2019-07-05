'''
(c) University of Liverpool 2019

All rights reserved.

@author: neilswainston
'''
# pylint: disable=too-many-public-methods
import unittest

from subliminal import math_utils


class Test(unittest.TestCase):
    '''Test class for math_utils.'''

    def test_isclose(self):
        '''Tests isclose method.'''
        self.assertTrue(math_utils.isclose(18.015054684, 18.015054684))
        self.assertFalse(math_utils.isclose(18.015054684, 18.015054))

    def test_linprog(self):
        '''Tests linprog method.'''
        a_matrix = [[-238, -142, 1, 248, 131, -1, -2, 1, 2],
                    [-144, -86, 0, 150, 80, 0, 0, 0, 0],
                    [-2, 0, 0, 2, 0, 0, 0, 0, 0],
                    [-57, -9, 0, 62, 4, 0, -1, 0, 1],
                    [-2, -1, 0, 2, 1, 0, 0, 0, 0],
                    [2.0, 1.0, 1.0, -2.0, -2.0, -1, 0, 1, 0]]

        bounds = [(1, 8), (1, 8), (1, 8), (1, 8), (1, 8), (0, 8), (0, 8),
                  (0, 8), (0, 8)]

        success, primals = math_utils.linprog([1] * len(a_matrix[0]),
                                              a_matrix,
                                              [0] * len(a_matrix),
                                              bounds)
        self.assertTrue(success)
        self.assertEqual(primals, [1, 1, 1, 1, 1, 0, 0, 0, 0])


if __name__ == "__main__":
    # import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
