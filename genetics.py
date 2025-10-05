from bayes import BayesianNetwork
from factor import *


class FamilyMember:
    """A single member of a family tree."""

    def __init__(self, name, sex, mother, father):
        """
        Parameters
        ----------
        name : str
            The name of the family member.
        sex : str
            The sex of the family member ("male" or "female")
        mother : FamilyMember
            The mother of the family member (or None if unknown)
        father : FamilyMember
            The father of the family member (or None if unknown)
        """

        self.name = name
        self.sex = sex
        self.mother = mother
        self.father = father

    def get_name(self):
        """Returns the name of the family member."""
        return self.name

    def get_sex(self):
        """Returns the sex of the family member."""
        return self.sex


class Male(FamilyMember):
    """A male family member."""

    def __init__(self, name, mother=None, father=None):
        super().__init__(name, "male", mother, father)


class Female(FamilyMember):
    """A female family member."""

    def __init__(self, name, mother=None, father=None):
        super().__init__(name, "female", mother, father)


def romanoffs():
    """A simple example of a family, using four members of the Russian royal family (the Romanoffs)."""
    alexandra = Female("alexandra")
    nicholas = Male("nicholas")
    alexey = Male("alexey", mother=alexandra, father=nicholas)
    anastasia = Female("anastasia", mother=alexandra, father=nicholas)
    return alexandra, nicholas, alexey, anastasia


def create_variable_domains(family):
    """Creates a dictionary mapping each variable to its domain, for the hemophilia network.

    For each family member, we create either 3 or 4 variables (3 if they’re male, 4 if they’re female).
    If N is the name of the family member, then we create the following variables:
        M_N: N’s maternally inherited gene
        P_N: N’s paternally inherited gene (if N is female)
        G_N: the genotype of N
        H_N: whether N has hemophilia

    The variables should be mapped to the following domains:
        - M_N: ['x', 'X']
        - P_N: ['x', 'X']
        - G_N: ['xx', 'xX', 'XX']
        - H_N: ['-', '+']

    Parameters
    ----------
    family : list[FamilyMember]
        the list of family members

    Returns
    -------
    dict[str, list[str]]
        a dictionary mapping each variable to its domain (i.e. its possible values)
    """
    # TODO: Implement this for Question Six.
    """pseudocode:
    begin w/empty dict
    iterate through the family members in family
        for each member of the family
            regardless of sex, they should have the maternally inherited gene
                add M_family_member to the dict
                map M_family_member to ['x', 'X']
            check if they're male or female
                if male:
                    overall here, I think the guidelines are wrong? 
                    The male genotype should include a y chromosome, which impacts the probability of getting hemophilia
                    add G_family_member and H_family_member to the dict
                    map G_family_member to ['xy', 'Xy']
                    map H_family_member to ['-', '+']
                if female:
                    add P_family_member, G_family_member, and H_family_member to the dict
                    map P_family_member to ['x', 'X']
                    map G_family_member to ['xx', 'xX', 'XX']
                    map H_family_member to ['-', '+']
    return the dict
    """
    variable_domains = {}
    for member in family:
        name = member.get_name() # uses the get_name() method in FamilyMember class
        variable_domains[f"M_{name}"] = ['x', 'X'] # you always have an X chromosome (maternally inherited)
        if member.get_sex() == "male":
            variable_domains[f"G_{name}"] = ['xy', 'Xy'] # males have one X and one Y chromosome
            variable_domains[f"H_{name}"] = ['-', '+']
        elif member.get_sex() == "female":
            variable_domains[f"P_{name}"] = ['x', 'X']
            variable_domains[f"G_{name}"] = ['xx', 'xX', 'XX']
            variable_domains[f"H_{name}"] = ['-', '+']

    return variable_domains


def create_hemophilia_cpt(person):
    """Creates a conditional probability table (CPT) specifying the probability of hemophilia, given one's genotype.

    Parameters
    ----------
    person : FamilyMember
        the family member whom the CPT pertains to

    Returns
    -------
    Factor
        a Factor specifying the probability of hemophilia, given one's genotype
    """
    # TODO: Implement this for Question Seven.
    """pseudocode:

    this is pretty brute-force. there's probably a way we can do fun itertools stuff, but I'm having trouble
    with the fact that create_variable_domains() takes in a family, not a person.


    a Factor needs a variable list and a values dict, with our variables being G_person and H_person
    initialize variable list as ["G_person", "H_person"]
    initialize values dict as empty dict

    person only gives us a name and a sex

    check the sex of the person
    if the person's a male, the possible genotypes are "xy" and "Xy"
        add to values_dict:
            ("xy", "+") -> 0.0
            ("xy", "-") -> 1.0
            ("Xy", "+") -> 1.0
            ("Xy", "-") -> 0.0
    if the person's a female, the possible genotypes are "xx", "xX", and "XX"
        add to values_dict:
            ("xx", "+") -> 0.0
            ("xx", "-") -> 1.0

            ("xX", "+") -> 0.0
            ("xX", "-") -> 1.0

            ("XX", "-") -> 0.0
            ("XX", "+") -> 1.0

    return a Factor with the variable list and values dict
    """
    name = person.get_name()
    variable_list = [f"G_{name}", f"H_{name}"]
    values_dict = {}

    # a male needs only one uppercase X chromosome to have hemophilia
    if person.get_sex() == "male":
        values_dict[("xy", "+")] = 0.0
        values_dict[("xy", "-")] = 1.0

        values_dict[("Xy", "+")] = 1.0
        values_dict[("Xy", "-")] = 0.0

    # a female needs two uppercase X chromosomes to have hemophilia
    elif person.get_sex() == "female":
        values_dict[("xx", "+")] = 0.0
        values_dict[("xx", "-")] = 1.0

        values_dict[("xX", "+")] = 0.0
        values_dict[("xX", "-")] = 1.0

        values_dict[("XX", "+")] = 1.0
        values_dict[("XX", "-")] = 0.0

    return Factor(variable_list, values_dict)


def create_genotype_cpt(person):
    """Creates a conditional probability table (CPT) specifying the probability of a genotype, given one's inherited genes.

    Parameters
    ----------
    person : FamilyMember
        the family member whom the CPT pertains to

    Returns
    -------
    Factor
        a Factor specifying the probability of a genotype, given one's inherited genes
    """
    # TODO: Implement this for Question Eight.
    """
    pseudocode:
    check the sex of the person
        if male, variable list is ["M_person", "G_person"]
        if female, variable list is ["M_person", "P_person", "G_person"]

    initialize empty values dict

    if male:
        possible maternally inherited genes: "x", "X"
        possible genotypes: "xy", "Xy"
        add to values_dict:
            ("x", "xy") -> 1.0
            ("x", "Xy") -> 0.0
            ("X", "xy") -> 0.0
            ("X", "Xy") -> 1.0
    elif female:
        possible maternally inherited genes: "x", "X"
        possible paternally inherited genes: "x", "X"
        possible genotypes: "xx", "xX", "XX"
        add to values_dict:
            ("x", "x", "xx") -> 1.0
            ("x", "x", "xX") -> 0.0
            ("x", "x", "XX") -> 0.0

            ("x", "X", "xx") -> 0.0
            ("x", "X", "xX") -> 1.0
            ("x", "X", "XX") -> 0.0

            ("X", "x", "xx") -> 0.0
            ("X", "x", "xX") -> 1.0
            ("X", "x", "XX") -> 0.0

            ("X", "X", "xx") -> 0.0
            ("X", "X", "xX") -> 0.0
            ("X", "X", "XX") -> 1.0

    return a Factor with the variable list and values dict
    """

    name = person.get_name()
    values_dict = {}

    if person.get_sex() == "male":
        variable_list = [f"M_{name}", f"G_{name}"]
        values_dict[("x", "xy")] = 1.0
        values_dict[("x", "Xy")] = 0.0

        values_dict[("X", "xy")] = 0.0
        values_dict[("X", "Xy")] = 1.0

    elif person.get_sex() == "female":
        variable_list = [f"M_{name}", f"P_{name}", f"G_{name}"]
        values_dict[("x", "x", "xx")] = 1.0
        values_dict[("x", "x", "xX")] = 0.0
        values_dict[("x", "x", "XX")] = 0.0

        values_dict[("x", "X", "xx")] = 0.0
        values_dict[("x", "X", "xX")] = 1.0
        values_dict[("x", "X", "XX")] = 0.0

        values_dict[("X", "x", "xx")] = 0.0
        values_dict[("X", "x", "xX")] = 1.0
        values_dict[("X", "x", "XX")] = 0.0

        values_dict[("X", "X", "xx")] = 0.0
        values_dict[("X", "X", "xX")] = 0.0
        values_dict[("X", "X", "XX")] = 1.0
    
    return Factor(variable_list, values_dict)

def create_maternal_inheritance_cpt(person):
    """Creates a conditional probability table (CPT) specifying the probability of the gene inherited from one's mother.

    Parameters
    ----------
    person : FamilyMember
        the family member whom the CPT pertains to

    Returns
    -------
    Factor
        a Factor specifying the probability of the gene inherited from the family member's mother.
    """
    # TODO: Implement this for Question Nine.
    """
    pseudocode:
    we've got two kinds of people here:
        those with no known mother
            in this case, the variable list is just ["M_person"]
            values dict is:
                ("x") -> 29999/30000
                ("X") -> 1/30000
        those with a known mother
            in this case, the variable list is ["G_mother", "M_person"]
            values dict is:
                ("xx", "x") -> 1.0
                ("xx", "X") -> 0.0

                ("xX", "x") -> 0.5
                ("xX", "X") -> 0.5

                ("XX", "x") -> 0.0
                ("XX", "X") -> 1.0
    return a Factor with the variable list and values dict
    """
    name = person.get_name()
    values_dict = {}
    if person.mother is None:
        variable_list = [f"M_{name}"]
        values_dict[("x",)] = 29999/30000
        values_dict[("X",)] = 1/30000
    else:
        mother_name = person.mother.get_name()
        variable_list = [f"G_{mother_name}", f"M_{name}"]
        values_dict[("xx", "x")] = 1.0
        values_dict[("xx", "X")] = 0.0

        values_dict[("xX", "x")] = 0.5
        values_dict[("xX", "X")] = 0.5

        values_dict[("XX", "x")] = 0.0
        values_dict[("XX", "X")] = 1.0

    return Factor(variable_list, values_dict)

def create_paternal_inheritance_cpt(person):
    """Creates a conditional probability table (CPT) specifying the probability of the gene inherited from one's father.

    Parameters
    ----------
    person : FamilyMember
        the family member whom the CPT pertains to

    Returns
    -------
    Factor
        a Factor specifying the probability of the gene inherited from the family member's father.
    """
    # TODO: Implement this for Question Ten.
    """
    pseudocode:
    three kinds of people here:
        men:
            variable list is ["P_person"]
            values dict is:
                ("y") -> 1.0
        women with no known father:
            variable list is ["P_person"]
            values dict is:
                ("x") -> 29999/30000
                ("X") -> 1/30000
        women with a known father:
            variable list is ["G_father", "P_person"]
            values dict is:
                ("xy", "x") -> 1.0
                ("xy", "X") -> 0.0

                ("Xy", "x") -> 0.0
                ("Xy", "X") -> 1.0
    return a Factor with the variable list and values dict
    """
    name = person.get_name()
    values_dict = {}
    if person.get_sex() == "male":
        variable_list = [f"P_{name}"]
        values_dict[("y",)] = 1.0
            # need to distinguish between the string "x" and the single-value tuple ("x",), as
            # values_dict keys need to be tuples (see factor.py, line 15)
    elif person.get_sex() == "female":
        if person.father is None:
            variable_list = [f"P_{name}"]
            values_dict[("x",)] = 29999/30000 
            values_dict[("X",)] = 1/30000
        else:
            father_name = person.father.get_name()
            variable_list = [f"G_{father_name}", f"P_{name}"]
            values_dict[("xy", "x")] = 1.0
            values_dict[("xy", "X")] = 0.0

            values_dict[("Xy", "x")] = 0.0
            values_dict[("Xy", "X")] = 1.0

    return Factor(variable_list, values_dict)

def create_family_bayes_net(family):
    """Creates a Bayesian network that models the genetic inheritance of hemophilia within a family.

    Parameters
    ----------
    family : list[FamilyMember]
        the members of the family

    Returns
    -------
    BayesianNetwork
        a Bayesian network that models the genetic inheritance of hemophilia within the specified family
    """
    domains = create_variable_domains(family)
    cpts = []
    for person in family:
        if person.get_sex() == "female":
            cpts.append(create_paternal_inheritance_cpt(person))
        cpts.append(create_maternal_inheritance_cpt(person))
        cpts.append(create_genotype_cpt(person))
        cpts.append(create_hemophilia_cpt(person))
    return BayesianNetwork(cpts, domains)
