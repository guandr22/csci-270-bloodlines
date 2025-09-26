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
    go through the entries in factor.values
    start copying over those entries to the new dictionary:
        if an entry shares the non-removed variable values with an existing entry,
            add the entry's value to the existing entry
        otherwise,
            create a new entry with the same variables, minus the variable we're removing
    remove the variable from factor.variables.
    """
    new_values = {}
    new_variables = [var for var in factor.variables if var != variable]

    # iterate through the entries in factor.values
    for event, value in factor.values.items():

        # create new event without the marginalized variable
        old_event_with_variables = zip(factor.variables, event)
        new_event = tuple(val for var, val in old_event_with_variables if var != variable) 
        if new_event in new_values:
            new_values[new_event] += value # consolidate the probabilities
        else:
            new_values[new_event] = value

    return Factor(new_variables, new_values)


def multiply_factors(factors, domains):
    """TODO: Implement this for Question Three."""