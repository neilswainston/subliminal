'''
(c) University of Liverpool 2019

All rights reserved.

@author: neilswainston
'''
from collections import Counter
import itertools

from cobra.io import read_sbml_model, write_sbml_model

from subliminal import math_utils
from synbiochem.utils import chem_utils


def balance_sbml_model(in_filename, out_filename):
    '''Attempts to mass / charge balance SBML L3 fbc models.'''
    model = balance_model(read_sbml_model(in_filename))
    write_sbml_model(model, out_filename)


def balance_model(model, verbose=True):
    '''Attempts to mass / charge balance cobra models.'''
    for react in model.reactions:
        # If NOT exchange reaction:
        if len(react.metabolites) > 1 and len(check_react(react)) > 1:
            reac_def = [(met.formula, met.charge, stoich, met.id)
                        for met, stoich in react.metabolites.items()]

            # Get reaction compartment:
            cmpt = Counter([met.compartment
                            for met in react.metabolites]).most_common(1)[0][0]

            balanced, unchanged, new_reac_def = \
                balance_reac(reac_def, max_stoich=10)

            # If reaction has been 'fixed' update stoichiometries:
            if balanced and not unchanged:
                if verbose:
                    print('Reaction %s fixed:' % react.id)
                    print('FROM:\t%s' % react.build_reaction_string())

                # Remove existing metabolites:
                mets = {met.id: 0
                        for met, stoich
                        in react.metabolites.items()}

                react.add_metabolites(mets, combine=False)

                # Add updated metabolites:
                for val in new_reac_def:
                    try:
                        metabolite = [met for met in mets if met == val[3]]
                    except:
                        print(mets)

                    if metabolite:
                        mets[metabolite[0]] = val[2]
                    else:
                        # This makes 2 assumptions:
                        #
                        # 1) That h_c, h_m, h2o_c etc are already in the model.
                        # 2) That the correct compartment for H or H20 is
                        #    that of the most frequent compartment for the
                        #    other metabolites in the reaction.
                        mets[val[3] + '_' + cmpt] = val[2]

                react.add_metabolites(mets)

                if verbose:
                    print('TO:\t%s' % react.build_reaction_string())

            elif not balanced and verbose:
                print('Reaction %s unbalanced' % react.id)

    return model


def balance_reac(react_def, optional_comp=None, max_stoich=8.0):
    '''Applies linear programming to balance reaction.'''
    if optional_comp is None:
        # Formula, charge, metabolite id (assume BiGG id)
        optional_comp = [('H', 1, 'h'),
                         ('H2O', 0, 'h2o')]

    formulae = [[x[0] for x in react_def if x[2] <= 0],
                [x[0] for x in react_def if x[2] > 0]]
    formulae.extend([[x[0] for x in optional_comp] for _ in range(2)])

    charges = [[x[1] for x in react_def if x[2] <= 0],
               [x[1] for x in react_def if x[2] > 0]]
    charges.extend([[x[1] for x in optional_comp] for _ in range(2)])

    ids = [[x[3] for x in react_def if x[2] <= 0],
           [x[3] for x in react_def if x[2] > 0]]
    ids.extend([[x[2] for x in optional_comp] for _ in range(2)])

    all_elem_comp = [[_get_elem_comp(formula, idx) for formula in formulae]
                     for idx, formulae in enumerate(formulae)]

    (balanced, stoichs) = _optimise(_get_elem_matrix(all_elem_comp,
                                                     charges),
                                    max_stoich, len(optional_comp))

    balanced_def = sorted(_get_react_def(stoichs, formulae, charges, ids)) \
        if balanced else react_def

    return balanced, \
        balanced and _compare_react_defs(react_def, balanced_def), \
        balanced_def


def check_model(model):
    '''Checks a cobra model for mass / charge balancing.'''
    model_balance = {}

    for reaction in model.reactions:
        # If NOT exchange reaction:
        if len(reaction.metabolites) > 1:
            # Check mass balance and remove tiny coefficients due to rounding:
            reaction_balance = check_react(reaction)

            if reaction_balance:
                model_balance[reaction.id] = reaction_balance

    return model_balance


def check_react(reaction):
    '''Checks a cobra reaction for mass / charge balancing.'''
    if len(reaction.metabolites) > 1:
        # Check mass balance and remove tiny coefficients due to rounding:
        mass_balance = reaction.check_mass_balance()

        return {key: val for key, val in mass_balance.items()
                if abs(val) > 1e-6}

    raise ValueError('No reaction metabolites')


def _get_react_def(stoichs, all_formulae, all_charges, all_ids):
    '''Formats the input into (formula, charge, stoichiometry) tuples.'''
    return [(a[0], b, a[1] * c, d)
            for a, b, c, d in zip([(x, -1 if idx % 2 == 0 else 1)
                                   for idx, formulae in enumerate(all_formulae)
                                   for x in formulae],
                                  [y for charges in all_charges
                                   for y in charges],
                                  _simplify_stoichs(stoichs),
                                  [z for ids in all_ids for z in ids])
            if c > 1e-8]


def _compare_react_defs(def1, def2):
    '''Returns True/False depending on whether reaction definitions are
    equal.'''
    return def1 is not None and def2 is not None and \
        len(def1) == len(def2) and \
        all([x[0] == y[0] and x[1] == y[1] and math_utils.isclose(x[2], y[2])
             for x, y in zip(sorted(def1), sorted(def2))])


def _simplify_stoichs(stoichs):
    '''Attempts to simplify stoichs of 1.00000001 to 1.0.'''
    return [float(round(stoich)) if math_utils.isclose(stoich, round(stoich))
            else stoich for stoich in stoichs]


def _get_elem_comp(formula, idx):
    '''Returns elemental compositions, multiplying stoichiometry by -1
    if reactants.'''
    elem_comp = chem_utils.get_elem_comp(formula)
    elem_comp.update((x, y * (-1 if idx % 2 == 0 else 1))
                     for x, y in elem_comp.items())
    return elem_comp


def _get_elem_matrix(all_elem_comp, all_charges):
    '''Gets the elemental (and charge) matrix representing rows of elements
    and columns of compounds.'''
    a_matrix = []
    elements = [elem_comp.keys() for elem_comps in all_elem_comp
                for elem_comp in elem_comps]

    for element in {item for sublist in elements for item in sublist}:
        a_matrix.append([elem_comp[element] if element in elem_comp else 0
                         for elem_comps in all_elem_comp
                         for elem_comp in elem_comps])

    a_matrix.append(list(itertools.chain(*([charge * (-1 if idx % 2 == 0
                                                      else 1)
                                            for charge in charges]
                                           for idx, charges
                                           in enumerate(all_charges)))))

    return a_matrix


def _optimise(a_matrix, max_stoich, num_opt_formulae):
    '''Optimise linear program and return result.'''
    bounds = [(1.0, max_stoich)] * (len(a_matrix[0]) -
                                    num_opt_formulae * 2) + \
        [(0.0, max_stoich)] * num_opt_formulae * 2

    return math_utils.linprog([1] * len(a_matrix[0]),
                              a_matrix,
                              [0.0] * (len(a_matrix)),
                              bounds)
