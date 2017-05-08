'''
synbiochem (c) University of Manchester 2016

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
import os
import unittest

import cobra
from subliminal import balance


class Test(unittest.TestCase):
    '''Test class.'''

    def test_balance(self):
        '''Tests balance method.'''
        filename = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                '../../models/MG1655_limonene.xml')
        model = balance.balance_model(cobra.io.read_sbml_model(filename))
        model_balance = balance.check_model_balance(model)
        self.assertEqual(len(model_balance), 1)

    def test_balance_unbalanced(self):
        '''Tests get_elem_comp method for unbalanced reaction.'''
        unbalanced = [('CO2', 0, -1.0, 'CO2'),
                      ('C5H7O4', -1, -1.0, 'C5H7O4'),
                      ('C3H3O3', -1, 1.0, 'C3H3O3')]

        is_balanced, was_balanced, balanced_def = \
            balance.balance_reac(unbalanced)

        self.assertTrue(is_balanced)
        self.assertFalse(was_balanced)
        self.assertEqual(sorted([('C3H3O3', -1, 2.0, 'C3H3O3'),
                                 ('C5H7O4', -1, -1.0, 'C5H7O4'),
                                 ('CO2', 0, -1.0, 'CO2'),
                                 ('H', 1, 1.0, 'CHEBI:24636')]),
                         sorted(balanced_def))

    def test_balance_balanced(self):
        '''Tests get_elem_comp method for balanced reaction.'''
        balanced = [('C5H7O4', -1, -1.0, 'C5H7O4'),
                    ('H', 1, 1.0, 'H'),
                    ('C3H3O3', -1, 2.0, 'C3H3O3'),
                    ('CO2', 0, -1.0, 'CO2')]

        is_balanced, was_balanced, balanced_def = \
            balance.balance_reac(balanced)

        self.assertTrue(is_balanced)
        self.assertTrue(was_balanced)
        self.assertEqual(sorted(balanced), sorted(balanced_def))

    def test_balance_problem(self):
        '''Tests get_elem_comp method for problematic reaction.'''
        problem = [('C144H238N2O57P2', -2.0, -1.0, 'C144H238N2O57P2'),
                   ('C86H142O9P', -1.0, -1.0, 'C86H142O9P'),
                   ('H', 1.0, 1.0, 'H'),
                   ('C150H248N2O62P2', -2.0, 1.0, 'C150H248N2O62P2'),
                   ('C80H131O4P', -2.0, 1.0, 'C80H131O4P')]

        is_balanced, was_balanced, balanced_def = \
            balance.balance_reac(problem, optional_comp=[])

        self.assertTrue(is_balanced)
        self.assertTrue(was_balanced)
        self.assertEqual(sorted(problem), sorted(balanced_def))


if __name__ == "__main__":
    unittest.main()
