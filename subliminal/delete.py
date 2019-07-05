'''
(c) University of Liverpool 2019

All rights reserved.

@author: neilswainston
'''
import sys

import cobra
from cobra.flux_analysis import single_gene_deletion

import subliminal.solve as solve


def delete(model, objectives):
    '''Perform single deletions.'''
    solve.set_objective(model, objectives, min_bound=1e-3)
    solve.solve(model)
    max_flux = model.solution.f
    return max_flux, single_gene_deletion(model)


def main(args):
    '''main method.'''
    val = iter(args[1:])
    objectives = {reac_id: float(obj_coeff)
                  for reac_id, obj_coeff in zip(val, val)}
    max_flux, result = delete(cobra.io.read_sbml_model(args[0]), objectives)

    for gene, flux in result[0].iteritems():
        percent_flux = (flux if flux > 1e-6 else 0.0) / max_flux
        print(gene, '{0:.4f}'.format(percent_flux))


if __name__ == '__main__':
    main(sys.argv[1:])
