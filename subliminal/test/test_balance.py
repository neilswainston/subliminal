'''
(c) University of Liverpool 2019

All rights reserved.

@author: neilswainston
'''
import os
import unittest

import cobra

from subliminal import balance


class Test(unittest.TestCase):
    '''Test class.'''

    def test_cobra_balance(self):
        '''Tests balance method.'''
        filename = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                '../../models/MG1655_limonene.xml')
        model = balance.balance_model(cobra.io.read_sbml_model(filename))
        model_balance = balance.check_model(model)
        self.assertEqual(len(model_balance), 2)

    def test_balance_unbalanced(self):
        '''Tests balance_reac method for unbalanced reaction.'''
        definition = [('CO2', 0, -1.0, 'CO2'),
                      ('C5H7O4', -1, -1.0, 'C5H7O4'),
                      ('C3H3O3', -1, 1.0, 'C3H3O3')]

        is_balanced, was_balanced, balanced_def = \
            balance.balance_reac(definition)

        self.assertTrue(is_balanced)
        self.assertFalse(was_balanced)
        self.assertEqual(sorted([('C3H3O3', -1, 2.0, 'C3H3O3'),
                                 ('C5H7O4', -1, -1.0, 'C5H7O4'),
                                 ('CO2', 0, -1.0, 'CO2'),
                                 ('H', 1, 1.0, 'h')]),
                         sorted(balanced_def))

    def test_balance_balanced(self):
        '''Tests balance_reac method for balanced reaction.'''
        definition = [('C5H7O4', -1, -1.0, 'C5H7O4'),
                      ('H', 1, 1.0, 'H'),
                      ('C3H3O3', -1, 2.0, 'C3H3O3'),
                      ('CO2', 0, -1.0, 'CO2')]

        is_balanced, was_balanced, balanced_def = \
            balance.balance_reac(definition)

        self.assertTrue(is_balanced)
        self.assertTrue(was_balanced)
        self.assertEqual(sorted(definition), sorted(balanced_def))

    def test_balance_problem(self):
        '''Tests balance_reac method for problematic reaction.'''
        definition = [('C144H238N2O57P2', -2.0, -1.0, 'C144H238N2O57P2'),
                      ('C86H142O9P', -1.0, -1.0, 'C86H142O9P'),
                      ('H', 1.0, 1.0, 'H'),
                      ('C150H248N2O62P2', -2.0, 1.0, 'C150H248N2O62P2'),
                      ('C80H131O4P', -2.0, 1.0, 'C80H131O4P')]

        is_balanced, was_balanced, balanced_def = \
            balance.balance_reac(definition, optional_comp=[])

        self.assertTrue(is_balanced)
        self.assertTrue(was_balanced)
        self.assertEqual(sorted(definition), sorted(balanced_def))

    def test_balance_unfeasible(self):
        '''Tests balance_reac method for unfeasible reaction.'''
        definition = [['C2H3O2', -1, 1.0, 'ac'],
                      ['H', 1, 1.0, 'h'],
                      ['HS2O3', -2, -1.0, 'tsul'],
                      ['C5H9NO4', 0, -1.0, 'acser'],
                      ['C3H6NO5S2', 0, 1.0, 'scys__L']]

        is_balanced, was_balanced, balanced_def = \
            balance.balance_reac(definition, optional_comp=[])

        self.assertFalse(is_balanced)
        self.assertFalse(was_balanced)
        self.assertEqual(sorted(definition), sorted(balanced_def))


if __name__ == '__main__':
    unittest.main()
