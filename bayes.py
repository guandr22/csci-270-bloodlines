from factor import *
from util import compute_elimination_order

class BayesianNetwork:
    """Represents a Bayesian network by its factors, i.e. the conditional probability tables (CPTs).

    Parameters
    ----------
    factors : list[factor.Factor]
        The factors of the Bayesian network
    domains : dict[str, list[str]]
        A dictionary mapping each variable to its possible values
    """

    def __init__(self, factors, domains):
        self.factors = factors
        self.domains = domains
        self.variables = set()
        for factor in self.factors:
            self.variables = self.variables | set(factor.variables)

    def __str__(self):
        return "\n\n".join([str(factor) for factor in self.factors])


def eliminate(bnet, variable):
    """Eliminates a variable from the Bayesian network.

    By "eliminate", we mean that the factors containing the variable are multiplied,
    and then the variable is marginalized (summed) out of the resulting factor.

    Parameters
    ----------
    variable : str
        the variable to eliminate from the Bayesian network

    Returns
    -------
    BayesianNetwork
        a new BayesianNetwork, equivalent to the current Bayesian network, after
        eliminating the specified variable
    """
    # TODO: Implement this for Question Four.

    # Find all factors with our target variable
    factors_with_variable = [factor for factor in bnet.factors if variable in factor.variables]

    # Multiply those factors together
    multiplied_factor = multiply_factors(factors_with_variable, bnet.domains)

    # Marginalize out the target variable from the multiplied factor
    marginalized_factor = marginalize(multiplied_factor, variable)

    # Create a new list of factors for the new Bayesian network, excluding the old factors with the target variable
    new_factors = [factor for factor in bnet.factors if variable not in factor.variables]
    new_factors.append(marginalized_factor)

    return BayesianNetwork(new_factors, bnet.domains)


def compute_marginal(bnet, vars):
    """Computes the marginal probability over the specified variables.

    This method uses variable elimination to compute the marginal distribution.

    Parameters
    ----------
    vars : set[str]
        the variables that we want to compute the marginal over
    """
    # TODO: Implement this for Question Five.

    # get the elimination order and moral graph (see util.py; we use the moral graph later)
    elim_order, _ = compute_elimination_order(bnet) 
    revised_elim_order = [var for var in elim_order if var not in vars]
    for var in revised_elim_order:
        bnet = eliminate(bnet, var)
    return multiply_factors(bnet.factors, bnet.domains)

def compute_conditional(bnet, event, evidence):
    """Computes the conditional probability of an event given the evidence event."""
    # TODO: Implement this for Question Five.
