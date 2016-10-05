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


def set_objective(model, objectives, min_bound=0, max_bound=999999):
    '''Sets the objective.'''
    for reaction in model.reactions:
        reaction.objective_coefficient = 0

    for reac_id, obj_coeff in objectives.iteritems():
        model.reactions.get_by_id(reac_id).objective_coefficient = obj_coeff
        set_bounds(model, {reac_id: [min_bound, max_bound]})


def solve(model, pfba=True):
    '''Solves the model, minimising total flux (if possible).'''
    if pfba and len(model.objective) == 1:
        cobra.flux_analysis.parsimonious.optimize_minimal_flux(model)
    else:
        model.optimize()


def print_solution(model):
    '''Print Solution.'''
    for reaction in model.reactions:
        if abs(reaction.x) > 1e-6:
            print reaction.id + '\t' + \
                reaction.build_reaction_string() + '\t' + \
                str(reaction.x)
