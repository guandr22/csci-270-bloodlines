from itertools import product
import itertools

class Factor:
    def __init__(self, variables, values):
        self.variables = variables
        self.values = values

    def __getitem__(self, event):
        key = []
        for var in self.variables:
            if var not in event:
                raise KeyError(f"Variable {var} not found in given event.")
            key.append(event[var])
        if tuple(key) in self.values:
            return self.values[tuple(key)]
        else:
            raise KeyError(f"No value assigned to event {event}.")

    def __str__(self):
        result = f"{self.variables}:"
        for event, value in self.values.items():
            result += f"\n  {event}: {value}"
        return result

    __repr__ = __str__


def events(vars, domains):
    """TODO: Implement this for Question One."""
    combinations = []
    for values in itertools.product(*(domains[var] for var in vars)):
        assignment = dict(zip(vars, values))
        combinations.append(assignment)
    return combinations
    
def marginalize(factor, variable):
    """TODO: Implement this for Question Two."""
    """
    Takes in a Factor object and a specified variable.
    Returns a new Factor object without that specified variable and, 
    accordingly, with all remaining joint probabilities consolidated
    in response to missing one variable.
    """
    
    """
    pseudocode: 

    start with a new empty dictionary 
    and a new list of variables without the variable we're removing
    go through the entries in factor.values
    start copying over those entries to the new dictionary:
        if an entry shares the non-removed variable values with an existing entry,
            add the entry's value to the existing entry
        otherwise,
            create a new entry with the same variables, minus the variable we're removing
    return a new factor object with the new list of variables and the new dictionary
    """
    new_values = {}
    new_variables = [var for var in factor.variables if var != variable]

    # iterate through the entries in factor.values
    for event, value in factor.values.items():

        # create new event without the marginalized variable
        old_event_with_variables = zip(factor.variables, event)

        # extracts only the values for the non-marginalized variables
        new_event = tuple(val for var, val in old_event_with_variables if var != variable) 
        if new_event in new_values:
            new_values[new_event] += value # consolidate the probabilities
        else:
            new_values[new_event] = value

    return Factor(new_variables, new_values)


def multiply_factors(factors, domains):
    """TODO: Implement this for Question Three."""
    """
    Takes in a list of two factors and a dictionary mapping variables to domains. 
    Produces a new factor wherein:
        1. For example, if variables x and y are in Factor 1 and variables y and z are in Factor 2,
        the new factor will have variables x, y, and z.
            We can generalize this to more than two variables in Factor 1 and Factor 2.
        2. The probability value of a given combination of x, y, and z will be the product
        of the probabilities of the combinations of x and y in Factor 1 and y and z in Factor 2.
            We can also generalize this to more than three variables in our new factor.
    """

    # get the list of all variables in the new factor
    new_variables = []
    for factor in factors:
        for var in factor.variables:
            if var not in new_variables:
                new_variables.append(var)

    # create a new dictionary to hold the new values
    new_values = {}

    # iterate through all possible events for the new factor
    for event in events(new_variables, domains):
        prob = 1
        for factor in factors:
            # extract the relevant part of the event for the current factor
            sub_event = {var: event[var] for var in factor.variables}
            prob *= factor[sub_event]
        new_values[tuple(event[var] for var in new_variables)] = prob 
        # add the new event and its probability to the new values dictionary
        # not going to lie, though, co-pilot did all the weird Python shit for me here - 
        # I only really know what's going on conceptually

    return Factor(new_variables, new_values)
