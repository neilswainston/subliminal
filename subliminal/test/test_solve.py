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

    def test_atp_solve(self):
        '''Tests solve method.'''
        model = _get_model()

        # Change glucose uptake:
        solve.set_bounds(model, {'EX_glc__D_e': [-1, None]})

        # Change the objective to ATPM and solve:
        solve.set_objective(model, 'ATPM')
        solve.solve(model)

        self.assertAlmostEqual(model.solution.f, 23.500, 3)

    def test_limonene_solve(self):
        '''Tests solve method.'''
        model = _get_model()

        # Change the objective to limonene excretion and solve:
        solve.set_objective(model, 'EX_MNXM956_c')
        solve.solve(model)

        solve.print_solution(model)

        self.assertAlmostEqual(model.solution.f, 3.201, 3)


def _get_model():
    '''Gets the model.'''
    filename = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                            '../../models/MG1655_limonene.xml')
    return cobra.io.read_sbml_model(filename)

if __name__ == "__main__":
    unittest.main()
