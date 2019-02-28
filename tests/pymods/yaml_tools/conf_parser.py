'''
    The token recognised by the configuration parser are declared here.
    There are two kind of token:
    - the parameters are just variables with a defautl value. They are used by
    constraints.
    - the constraints represent an actual check of the data. They consist in a
    function returning a boolean. The docstring of the function is used as the
    description for the testconf_explorer. There is not point in putting
    value_type, inherited, apply_to or use_params in the docstring but you
    should explain when the test fail and when it succeed.
'''
from __future__ import print_function, division, unicode_literals
from numpy import ndarray, isnan
from numpy.linalg import norm
from .meta_conf_parser import ConfParser
from .structures import Tensor

conf_parser = ConfParser()

# Parameters
conf_parser.parameter('tol_eq', default=1e-8, inherited=False)

tol_group = {
    'tol_rel', 'tol_abs', 'tol'
}


# Constraints
# default parameters for constraints are:
# value_type=float, inherited=True, apply_to='number' use_params=[], exclude={}
@conf_parser.constraint(exclude={'tol', 'ceil', 'ignore'})
def tol_rel(tol, ref, tested):
    '''
        Valid if the relative difference between the values is below the
        given tolerance.
    '''
    if abs(ref) + abs(tested) == 0.0:
        return True
    return abs(ref - tested) / (abs(ref) + abs(tested)) < tol


@conf_parser.constraint(exclude={'tol', 'ceil', 'ignore'})
def tol_abs(tol, ref, tested):
    '''
        Valid if the absolute difference between the values is below the
        given tolerance.
    '''
    return abs(ref - tested) < tol


@conf_parser.constraint(exclude={'ceil', 'tol_abs', 'tol_rel', 'ignore'})
def tol(tolv, ref, tested):
    '''
        Valid if both relative and absolute differences between the values
        are below the given tolerance.
    '''
    if abs(ref) + abs(tested) == 0.0:
        return True
    return (abs(ref - tested) / (abs(ref) + abs(tested)) < tolv and
            abs(ref - tested) < tolv)


@conf_parser.constraint(exclude={'tol', 'tol_abs', 'tol_rel', 'ignore'})
def ceil(ceil_val, ref, tested):
    '''
        Valid if the absolute value of the tested file tested is below the
        given tolerance.
    '''
    return abs(tested) < ceil_val


@conf_parser.constraint(exclude={'ceil', 'tol', 'tol_rel', 'tol_abs'})
def ignore(yes, ref, tested):
    '''
        Override numbers tests and always return the same result.
    '''
    return yes


@conf_parser.constraint(value_type=str, inherited=False, apply_to='this',
                        use_params=['tol_eq'])
def equation(eq, ref, tested, tol_eq):
    '''
        If the given expression return a number, its absolute value is compared
        to the optional tol_eq parameter. If the expression return a vector,
        it is the euclidian norm that is used.
    '''
    res = eval(eq, {}, {'this': tested, 'ref': ref})
    if isinstance(res, ndarray):
        return norm(res) < tol_eq
    else:
        return abs(res) < tol_eq


@conf_parser.constraint(value_type=list, inherited=False, apply_to='this',
                        use_params=['tol_eq'])
def equations(eqs, ref, tested, tol_eq):
    '''
        See equation. Same thing with a list of equations.
    '''
    for eq in eqs:
        res = eval(eq, {}, {'this': tested, 'ref': ref})
        if isinstance(res, ndarray):
            if norm(res) >= tol_eq:
                return False
        else:
            if abs(res) >= tol_eq:
                return False


@conf_parser.constraint(value_type=bool, inherited=True, apply_to='number')
def allow_nan(yes, ref, tested):
    '''
        If value is true, fail if a NaN value is encountered.
    '''
    return yes or not isnan(tested)


@conf_parser.constraint(value_type=bool, inherited=False, apply_to=Tensor)
def tensor_is_symetric(sym, ref, tested):
    '''
        If value is true:
            valid if tested tensor is symetric
        If value is false:
            valid if tested tensor is anti symetric
    '''
    if sym:
        return tested.is_symetric()
    else:
        return tested.is_anti_symetric()
