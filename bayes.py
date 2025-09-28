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

    # get the elimination order and moral graph (see util.py; we don't use the moral graph here)
    elim_order, _ = compute_elimination_order(bnet) 

    # revise the elimination order such that we only eliminate variables 
    # not included in our marginal probability calculation
    revised_elim_order = [var for var in elim_order if var not in vars]
    for var in revised_elim_order:
        bnet = eliminate(bnet, var)

    # multiply the remaining factors to get our marginal distribution
    return multiply_factors(bnet.factors, bnet.domains)

def extract_probability(factor, event):
    """Helper function that extracts the probability of a specific event from a factor.

    Parameters
    ----------
    factor : factor.Factor
        the factor from which to extract the probability
    event : dict[str, str]
        a dictionary mapping variable names to values, representing the event whose probability we want to extract

    Returns
    -------
    float
        the probability of the specified event according to the given factor
    """

    vars_in_event = set(event.keys())
    if set(factor.variables) != vars_in_event:
        raise ValueError("The variables in the event must match the variables in the factor.")

    # create a tuple of values corresponding to the order of variables in the factor
    event_tuple = tuple(event[var] for var in factor.variables)

    return factor.values.get(event_tuple, 0.0) #0.0 if event_tuple not in factor.values

def compute_conditional(bnet, event, evidence):
    """Computes the conditional probability of an event given the evidence event."""
    # TODO: Implement this for Question Five.
    """
    event: dict mapping variable names to values (e.g. {"X": "A", "Y": "AB", "Z": "B"} in our example about blood.)
    evidence: another dict mapping variable names to values (e.g. {"Y": "AB"} in our example about blood.), serving
        as the event(s) whose probability will become the denominator in our conditional probability calculation.

    returns: float, the conditional probability P(event | evidence)

    Our strategy relies on the definition of conditional probability:
        We'll calculate the joint probability of the event AND the evidence, then 
        divide that by the probability of the evidence.
    """

    # calculate the joint probability of event and evidence
    joint_event = {**event, **evidence} # merge the two dictionaries
    joint_vars = set(joint_event.keys())
    joint_marginal = compute_marginal(bnet, joint_vars)

    joint_prob = extract_probability(joint_marginal, joint_event)
    

    # calculate the marginal probability of evidence
    evidence_vars = set(evidence.keys())
    evidence_marginal = compute_marginal(bnet, evidence_vars)

    evidence_prob = extract_probability(evidence_marginal, evidence)

    if evidence_prob == 0.0:
        return 0.0  # avoid division by zero

    # apply the definition of conditional probability
    return joint_prob / evidence_prob