import numpy as np
import random

bases = ['A', 'C', 'G', 'T']
transitions = {
    'A': ['G'],
    'C': ['T'],
    'G': ['A'],
    'T': ['C']
}
transversions = {
    'A': ['C', 'T'],
    'C': ['A', 'G'],
    'G': ['C', 'T'],
    'T': ['A', 'G']
}

class MutationModel:
    def __init__(self, 
                 model_type: str,
                 base:str,
                 gamma: float = None,
                 alpha: float = None, 
                 beta: float = None, 
                 pi_a: float = None,
                 pi_c: float = None,
                 pi_g: float = None,
                 pi_t: float = None):
        """
        Initialize the mutation model.
        Parameters:
            model_type (str): Type of mutation model.
            base (str): The base to mutate from ('A', 'C', 'G', 'T').
            gamma (float): The mutation rate.
            alpha (float): Parameter for K2P and HKY85 models.
            beta (float): Parameter for K2P and HKY85 models.
            pi_a (float): Parameter for HKY85 and F81 models.
            pi_c (float): Parameter for HKY85 and F81 models.
            pi_g (float): Parameter for HKY85 and F81 models.
            pi_t (float): Parameter for HKY85 and F81 models.
        """
        
        self.model_type = model_type
        self.base = base
        self.alpha = alpha
        self.beta = beta
        self.gamma = gamma
        self.pi_a = pi_a
        self.pi_c = pi_c
        self.pi_g = pi_g
        self.pi_t = pi_t

        if self.gamma is None:
            raise ValueError("Gamma parameter must be provided for JC69 model.")
        if self.gamma < 0 or self.gamma > 1:
            raise ValueError("Gamma parameter must be between 0 and 1.")
        if self.base not in bases:
            raise ValueError(f"Base must be one of {bases}, got '{self.base}'.")
        if self.model_type not in ['random', 'JC69', 'K2P', 'F81', 'HKY85']:
            raise ValueError(f"Model type must be one of ['random', 'JC69', 'K2P', 'F81', 'HKY85'], got '{self.model_type}'.")
        
    def should_mutate(self):
        """
        Given the mutation rate, determine if a mutation should occur.
        Returns:
            bool: True if a mutation should occur, False otherwise.
        """
        return random.random() < self.gamma

        
    # define mutation models
    def random_mutation(self):
        """
        Randomly mutate a base to any of the other bases.
        """
        return random.choice(['A', 'C', 'G', 'T'])
    def JC69(self):
        """
        Jukes-Cantor model for equal base frequencies at a set mutation rate.
        """
        # Mutate to a different base
        possible_bases = [b for b in bases if b != self.base]
        return random.choice(possible_bases)
        
    def K2P(self):
        """
        Kimura 2-parameter model for transitions and transversions.
        Uses parameters alpha and beta for transition and transversion rates.
        """
        # Check if alpha and beta parameters are provided
        if self.alpha is None or self.beta is None:
            raise ValueError("Alpha and beta parameters must be provided for K2P model.")
        if self.alpha < 0 or self.beta < 0:
            raise ValueError("Alpha and beta parameters must be non-negative.")
        # Calculate the probabilities of transition and transversion
        transition_probability = self.alpha / (self.alpha + self.beta)
        transversion_probability = self.beta / (self.alpha + self.beta)
        # Mutate based on the probabilities
        if random.random() < transition_probability:
            # Transition mutation
            return random.choice(transitions[self.base])
        else:
            # Transversion mutation
            return random.choice(transversions[self.base])
        
    def F81(self):
        """
        Felsenstein 1981 model for nucleotide substitution.
        Uses parameters pi_a, pi_c, pi_g, pi_t for base frequencies.
        """
        # Check if base frequencies are provided
        if self.pi_a is None or self.pi_c is None or self.pi_g is None or self.pi_t is None:
            raise ValueError("Base frequencies (pi_a, pi_c, pi_g, pi_t) must be provided for F81 model.")
        if self.pi_a < 0 or self.pi_c < 0 or self.pi_g < 0 or self.pi_t < 0:
            raise ValueError("Base frequencies must be non-negative.")
        # Calculate the probabilities of mutation based on base frequencies
        total_freq = self.pi_a + self.pi_c + self.pi_g + self.pi_t
        if total_freq <= 0:
            raise ValueError("Base frequencies must sum to a positive value.")
        probabilities = {
            'A': self.pi_a / total_freq,
            'C': self.pi_c / total_freq,
            'G': self.pi_g / total_freq,
            'T': self.pi_t / total_freq
        }
        # Mutate based on the probabilities
        possible_bases = [b for b in bases if b != self.base]
        mutated_base = random.choices(possible_bases, weights=[probabilities[b] for b in possible_bases])[0]
        return mutated_base
        
    def HKY85(self):
        """
        Hasegawa-Kishino-Yano 1985 model for nucleotide substitution.
        Uses parameters pi_a, pi_c, pi_g, pi_t for base frequencies and alpha, beta for transition/transversion rates.
        """
        # Check if base frequencies and transition/transversion rates are provided
        if self.pi_a is None or self.pi_c is None or self.pi_g is None or self.pi_t is None:
            raise ValueError("Base frequencies (pi_a, pi_c, pi_g, pi_t) must be provided for HKY85 model.")
        if self.alpha is None or self.beta is None:
            raise ValueError("Alpha and beta parameters must be provided for HKY85 model.")
        if self.pi_a < 0 or self.pi_c < 0 or self.pi_g < 0 or self.pi_t < 0:
            raise ValueError("Base frequencies must be non-negative.")
        # Calculate the probabilities of mutation based on base frequencies
        total_freq = self.pi_a + self.pi_c + self.pi_g + self.pi_t
        if total_freq <= 0:
            raise ValueError("Base frequencies must sum to a positive value.")
        probabilities = {
            'A': self.pi_a / total_freq,
            'C': self.pi_c / total_freq,
            'G': self.pi_g / total_freq,
            'T': self.pi_t / total_freq
        }
        # Calculate the probabilities of transition and transversion
        transition_probability = self.alpha / (self.alpha + self.beta)
        transversion_probability = self.beta / (self.alpha + self.beta)

        # Mutate based on the probabilities
        if random.random() < transition_probability:
            # Transition mutation
            possible_bases = transitions[self.base]
            mutated_base = random.choices(possible_bases, weights=[probabilities[b] for b in possible_bases])[0]
        else:
            # Transversion mutation
            possible_bases = transversions[self.base]
            mutated_base = random.choices(possible_bases, weights=[probabilities[b] for b in possible_bases])[0]
        return mutated_base
    
    def mutate(self):
        """
        Simulate a mutation based on the specified model type.
        Returns the mutated base.
        """
        # Check if mutation should occur
        if not self.should_mutate():
            return self.base
        if self.model_type == 'random':
            return self.random_mutation()
        elif self.model_type == 'JC69':
            return self.JC69()
        elif self.model_type == 'K2P':
            return self.K2P()
        elif self.model_type == 'F81':
            return self.F81()
        elif self.model_type == 'HKY85':
            return self.HKY85()
        else:
            raise ValueError(f"Unknown model type: {self.model_type}")


from collections import Counter

def run_model_frequencies(model_name, model_kwargs, runs=10000):
    """
    Run a mutation model multiple times and collect mutation frequencies.
    
    Parameters:
        model_name (str): Model type.
        model_kwargs (dict): Keyword arguments for MutationModel constructor.
        runs (int): Number of simulations.
        
    Returns:
        dict: Frequencies of resulting mutations.
    """
    counts = Counter()
    for _ in range(runs):
        model = MutationModel(model_type=model_name, **model_kwargs)
        mutated = model.mutate()
        if mutated != model.base:
            counts[mutated] += 1
    total = sum(counts.values())
    return {k: v / total for k, v in counts.items()}


"""
Test the mutation models by running each one with a base and printing the frequencies of mutations.
def test_all_models():
    base = 'A'
    runs = 10000
    print(f"Testing each model with base '{base}' over {runs} runs\n")

    print("Random Model:")
    freqs = run_model_frequencies('random', {'base': base, 'gamma': 1.0}, runs)
    print(freqs)

    print("\nJC69 Model:")
    freqs = run_model_frequencies('JC69', {'base': base, 'gamma': 1.0}, runs)
    print(freqs)

    print("\nK2P Model:")
    freqs = run_model_frequencies('K2P', {
        'base': base, 'gamma': 1.0, 'alpha': 2.0, 'beta': 1.0
    }, runs)
    print(freqs)

    print("\nF81 Model:")
    freqs = run_model_frequencies('F81', {
        'base': base, 'gamma': 1.0,
        'pi_a': 0.1, 'pi_c': 0.2, 'pi_g': 0.3, 'pi_t': 0.4
    }, runs)
    print(freqs)

    print("\nHKY85 Model:")
    freqs = run_model_frequencies('HKY85', {
        'base': base, 'gamma': 1.0,
        'alpha': 2.0, 'beta': 1.0,
        'pi_a': 0.1, 'pi_c': 0.2, 'pi_g': 0.3, 'pi_t': 0.4
    }, runs)
    print(freqs)


# Run the test
test_all_models()
"""
