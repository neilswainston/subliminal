'''
(c) University of Liverpool 2019

All rights reserved.

@author: neilswainston
'''
import cobra


def set_bounds(model, bounds):
    '''Sets lower and upper bounds.'''
    for reac_id, bound_vals in bounds.items():
        if bound_vals[0] is not None:
            model.reactions.get_by_id(reac_id).lower_bound = bound_vals[0]
        if bound_vals[1] is not None:
            model.reactions.get_by_id(reac_id).upper_bound = bound_vals[1]


def set_objective(model, objectives, min_bound=0, max_bound=999999):
    '''Sets the objective.'''
    for reaction in model.reactions:
        reaction.objective_coefficient = 0

    for reac_id, obj_coeff in objectives.items():
        model.reactions.get_by_id(reac_id).objective_coefficient = obj_coeff
        set_bounds(model, {reac_id: [min_bound, max_bound]})


def solve(model, pfba=True):
    '''Solves the model, minimising total flux (if possible).'''
    if pfba:
        return cobra.flux_analysis.parsimonious.optimize_minimal_flux(model)

    return model.optimize()


def print_solution(model):
    '''Print Solution.'''
    for reaction in model.reactions:
        if abs(reaction.x) > 1e-6:
            print(reaction.id,
                  reaction.build_reaction_string(),
                  reaction.x)
