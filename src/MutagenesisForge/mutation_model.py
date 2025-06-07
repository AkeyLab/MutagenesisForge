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
                 gamma: float = 1.0,
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
            gamma (float): The mutation rate.
            alpha (float): Parameter for K2P and HKY85 models.
            beta (float): Parameter for K2P and HKY85 models.
            pi_a (float): Parameter for HKY85 and F81 models.
            pi_c (float): Parameter for HKY85 and F81 models.
            pi_g (float): Parameter for HKY85 and F81 models.
            pi_t (float): Parameter for HKY85 and F81 models.
        """
        
        self.model_type = model_type
        self.alpha = alpha
        self.beta = beta
        self.gamma = gamma
        self.pi_a = pi_a
        self.pi_c = pi_c
        self.pi_g = pi_g
        self.pi_t = pi_t

        if self.model_type not in ['random', 'JC69', 'K2P', 'F81', 'HKY85']:
            raise ValueError(f"Model type must be one of ['random', 'JC69', 'K2P', 'F81', 'HKY85'], got '{self.model_type}'.")
        
    def should_mutate(self) -> bool:
        """
        Given the mutation rate, determine if a mutation should occur.
        Returns:
            bool: True if a mutation should occur, False otherwise.
        """
        return random.random() < self.gamma

        
    # define mutation models
    def random_mutation(self) -> str:
        """
        Randomly mutate a base to any of the other bases.
        """
        return random.choice(['A', 'C', 'G', 'T'])
    def JC69(self, base) -> str:
        """
        Jukes-Cantor model for equal base frequencies at a set mutation rate.
        """
        # Mutate to a different base
        possible_bases = [b for b in bases if b != base]
        return random.choice(possible_bases)

    def K2P(self, base) -> str:
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
            return random.choice(transitions[base])
        else:
            # Transversion mutation
            return random.choice(transversions[base])

    def F81(self, base) -> str:
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
        possible_bases = [b for b in bases if b != base]
        mutated_base = random.choices(possible_bases, weights=[probabilities[b] for b in possible_bases])[0]
        return mutated_base
        
    def HKY85(self, base) -> str:
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
            possible_bases = transitions[base]
            mutated_base = random.choices(possible_bases, weights=[probabilities[b] for b in possible_bases])[0]
        else:
            # Transversion mutation
            possible_bases = transversions[base]
            mutated_base = random.choices(possible_bases, weights=[probabilities[b] for b in possible_bases])[0]
        return mutated_base
    
    def mutate(self, base:str) -> str:
        """
        Simulate a mutation based on the specified model type.
        Returns the mutated base.
        """
        if base not in bases:
            raise ValueError(f"Base must be one of {bases}, got '{base}'.")
        # Check if mutation should occur
        if not self.should_mutate():
            return base
        if self.model_type == 'random':
            return self.random_mutation()
        elif self.model_type == 'JC69':
            return self.JC69(base=base)
        elif self.model_type == 'K2P':
            return self.K2P(base=base)
        elif self.model_type == 'F81':
            return self.F81(base=base)
        elif self.model_type == 'HKY85':
            return self.HKY85(base=base)
        else:
            raise ValueError(f"Unknown model type: {self.model_type}")
    
    def get_mutation_probability(self, base:str, mutated_base:str) -> float:
        """
        Get the mutation probability for a given base to mutate into a specific base.
        """
        if base not in bases or mutated_base not in bases:
            raise ValueError(f"Both base and mutated_base must be one of {bases}.")

        if self.model_type == 'random':
            return self.gamma / 3
        elif self.model_type == 'JC69':
            return self.gamma / 3
        elif self.model_type == 'K2P':
            if base in transitions and mutated_base in transitions[base]:
                return self.gamma * self.alpha / (self.alpha + 2 * self.beta)
            elif base in transversions and mutated_base in transversions[base]:
                return self.gamma * self.beta / (self.alpha + 2 * self.beta)

        elif self.model_type == 'F81':
            total_freq = self.pi_a + self.pi_c + self.pi_g + self.pi_t
            if total_freq <= 0:
                raise ValueError("Base frequencies must sum to a positive value.")
            probabilities = {
                'A': self.pi_a / total_freq,
                'C': self.pi_c / total_freq,
                'G': self.pi_g / total_freq,
                'T': self.pi_t / total_freq
            }
            if mutated_base in probabilities:
                return self.gamma * probabilities[mutated_base]
        elif self.model_type == 'HKY85':
            total_freq = self.pi_a + self.pi_c + self.pi_g + self.pi_t
            if total_freq <= 0:
                raise ValueError("Base frequencies must sum to a positive value.")
            probabilities = {
                'A': self.pi_a / total_freq,
                'C': self.pi_c / total_freq,
                'G': self.pi_g / total_freq,
                'T': self.pi_t / total_freq
            }
            transition_probability = self.alpha / (self.alpha + self.beta)
            transversion_probability = self.beta / (self.alpha + self.beta)
            if base in transitions and mutated_base in transitions[base]:
                return self.gamma * transition_probability * probabilities[mutated_base]
            elif base in transversions and mutated_base in transversions[base]:
                return self.gamma * transversion_probability * probabilities[mutated_base]

        raise ValueError(f"Cannot calculate mutation probability for {base} to {mutated_base} with model {self.model_type}.")