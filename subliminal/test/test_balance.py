'''
synbiochem (c) University of Manchester 2016

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
import os
import unittest

import cobra

import subliminal.balance as balance


class Test(unittest.TestCase):
    '''Test class.'''

    def test_balance(self):
        '''Tests balance method.'''
        filename = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                '../../models/MG1655_limonene.xml')
        model = balance.balance(cobra.io.read_sbml_model(filename))
        model_balance = balance.check_model_balance(model)
        self.assertEqual(len(model_balance), 1)


if __name__ == "__main__":
    unittest.main()
