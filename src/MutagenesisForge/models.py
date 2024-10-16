import numpy as np

"""
this module contains evolutionary models for simulating mutations
"""

bases = ["A", "C", "G", "T"]


def random_mutation(base):
    """
    Given a base, return a random mutation

    Parameters:
        base (str): base to mutate

    Returns:
        str: mutated base
    """
    return np.random.choice([b for b in bases if b != base])


def K2P(base, alpha, beta):
    """
    Given a base, return a mutation based on the Kimura 2-parameter model probabilities

    Parameters:
        base (str): base to mutate
        alpha (float): transition probability
        beta (float): transversion probability

    Returns:
        str: mutated base

    Kimura 2-parameter model:

        A   C    G    T
    A    [0, alpha, beta, beta]
    C    [alpha, 0, beta, beta]
    G    [beta, beta, 0, alpha]
    T    [beta, beta, alpha, 0]

    """

    alpha = float(alpha)
    beta = float(beta)

    # ensure alpha + beta*2 = 1
    if not np.isclose(alpha + beta * 2, 1):
        raise ValueError("alpha + 2*beta must equal")

    # transition-transversion probabilities based on tstv ratio
    matrix_model = {
        "A": {"A": 0, "C": alpha, "G": beta, "T": beta},
        "C": {"A": alpha, "C": 0, "G": beta, "T": beta},
        "G": {"A": beta, "C": beta, "G": 0, "T": alpha},
        "T": {"A": beta, "C": beta, "G": alpha, "T": 0},
    }

    return np.random.choice(bases, p=[matrix_model[base][b] for b in bases])


def K3P(base, alpha, beta, gamma):
    """
    Given a base, return a mutation based on the Kimura 3-parameter model probabilities

    Parameters:
        base (str): base to mutate
        alpha (float): transition probability
        beta (float): transversion probability
        gamma (float): transversion probability

    Returns:
        str: mutated base

    Kimura 3-parameter model:

            A   C    G    T
        A    [0, alpha, beta, gamma]
        C    [alpha, 0, gamma, beta]
        G    [beta, gamma, 0, alpha]
        T    [gamma, beta, alpha, 0]
    """

    alpha = float(alpha)
    beta = float(beta)
    gamma = float(gamma)

    # ensure alpha + beta + gamma = 1
    if not np.isclose(alpha + beta + gamma, 1):
        raise ValueError("alpha + beta + gamma must equal 1")

    # transition-transversion probabilities
    matrix_model = {
        "A": {"A": 0, "C": alpha, "G": beta, "T": gamma},
        "C": {"A": alpha, "C": 0, "G": gamma, "T": beta},
        "G": {"A": beta, "C": gamma, "G": 0, "T": alpha},
        "T": {"A": gamma, "C": beta, "G": alpha, "T": 0},
    }

    return np.random.choice(bases, p=[matrix_model[base][b] for b in bases])
