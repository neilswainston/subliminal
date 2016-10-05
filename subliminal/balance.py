'''
synbiochem (c) University of Manchester 2016

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
from synbiochem.utils import chem_utils


def balance(model):
    '''Attempts to mass / charge balance models.'''
    for reaction in model.reactions:
        # If NOT exchange reaction:
        if len(reaction.metabolites) > 1 and \
                len(check_reaction_balance(reaction)) > 1:
            reac_def = [(met.formula, met.charge, stoich, met.id)
                        for met, stoich in reaction.metabolites.iteritems()]
            result = chem_utils.balance(reac_def, max_stoich=10)

            # If reaction has been 'fixed' update stoichiometries:
            if result[0] and not result[1]:
                print reaction.build_reaction_string()
                metabolites = {met: stoich
                               for met, stoich
                               in reaction.metabolites.iteritems()}
                reaction.clear_metabolites()

                for val in result[2]:
                    metabolite = [met for met in metabolites
                                  if met.id == val[3]][0]
                    metabolites[metabolite] = val[2]

                reaction.add_metabolites(metabolites)
                print reaction.build_reaction_string()
                print '\n'

    return model


def check_model_balance(model):
    '''Checks a model for mass / charge balancing.'''
    model_balance = {}

    for reaction in model.reactions:
        # If NOT exchange reaction:
        if len(reaction.metabolites) > 1:
            # Check mass balance and remove tiny coefficients due to rounding:
            reaction_balance = check_reaction_balance(reaction)

            if len(reaction_balance) > 0:
                model_balance[reaction.id] = reaction_balance

    return model_balance


def check_reaction_balance(reaction):
    '''Checks a reaction for mass / charge balancing.'''
    if len(reaction.metabolites) > 1:
        # Check mass balance and remove tiny coefficients due to rounding:
        mass_balance = reaction.check_mass_balance()

        return {key: val for key, val in mass_balance.iteritems()
                if val > 1e-6}
