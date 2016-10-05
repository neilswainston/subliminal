'''
synbiochem (c) University of Manchester 2016

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
import os
import unittest

import cobra

import subliminal.solve as solve


class Test(unittest.TestCase):
    '''Test class.'''

    def test_solve(self):
        '''Tests solve method.'''
        filename = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                '../../models/MG1655_limonene.xml')
        model = cobra.io.read_sbml_model(filename)

        # Remove duplicate reactions:
        model.remove_reactions(['MNXR341', 'MNXR55238', 'MNXR69289'])

        # Change glucose and oxygen uptake:
        solve.set_bounds(model, {'EX_glc__D_e': [-1, None],
                                 'EX_o2_e': [-float('inf'), None]})

        # Change the objective to ATPM and solve:
        solve.set_objective(model, 'ATPM')
        solve.solve(model)

        self.assertAlmostEqual(model.solution.f, 23.5)


if __name__ == "__main__":
    unittest.main()
