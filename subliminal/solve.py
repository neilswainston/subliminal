'''
synbiochem (c) University of Manchester 2016

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
import cobra


def set_bounds(model, bounds):
    '''Sets lower and upper bounds.'''
    for reac_id, bound_vals in bounds.iteritems():
        if bound_vals[0] is not None:
            model.reactions.get_by_id(reac_id).lower_bound = bound_vals[0]
        if bound_vals[1] is not None:
            model.reactions.get_by_id(reac_id).upper_bound = bound_vals[1]


def set_objective(model, reac_id):
    '''Sets the objective.'''
    model.objective = reac_id
    set_bounds(model, {reac_id: [0, 999999]})


def solve(model):
    '''Solves the model, minimising total flux.'''
    cobra.flux_analysis.parsimonious.optimize_minimal_flux(model)


def print_solution(model):
    '''Print Solution.'''
    for reaction in model.reactions:
        if abs(reaction.x) > 1e-6:
            print reaction.id + '\t' + \
                reaction.build_reaction_string() + '\t' + \
                str(reaction.x)
