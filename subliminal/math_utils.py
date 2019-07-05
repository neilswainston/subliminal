'''
(c) University of Liverpool 2019

All rights reserved.

@author: neilswainston
'''
import scipy.optimize


def isclose(value_1, value_2, rel_tol=1e-09, abs_tol=0.0):
    '''Compares floating point numbers.'''
    return abs(value_1 - value_2) <= max(rel_tol * max(abs(value_1),
                                                       abs(value_2)), abs_tol)


def linprog(c_vector, A_eq, b_eq, bounds):
    '''Solve linear programming problem with GLPK.'''
    res = scipy.optimize.linprog(c_vector, A_eq=A_eq, b_eq=b_eq, bounds=bounds)

    return res.status == 0, \
        [round(x, 8) for x in res.x] if res.status == 0 else None
