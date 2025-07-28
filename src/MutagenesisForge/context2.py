from contextlib import contextmanager
import gzip
import pysam
from Bio import SeqIO
import numpy as np
from collections import defaultdict
from functools import lru_cache
from typing import Tuple, List, Optional, Dict, Any
import logging
import os
import yaml
import random

from .mutation_model import MutationModel
from .utils import load_parameter_from_yaml, check_yaml_variable

"""
This module contains functions for creating a vcf file of random mutations.
"""

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Constants
VALID_BASES = {'A', 'C', 'G', 'T'}

# Codon to amino acid mapping (moved to module level for efficiency)
CODON_TO_AMINO = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M", "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T", "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*", "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K", "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W", "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R", "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G"
}


@contextmanager
def my_open(filename: str, mode: str):
    """A wrapper for open/gzip.open logic as a context manager"""
    try:
        with (gzip.open(filename, mode + "t") if filename.endswith(".gz") else open(filename, mode)) as open_file:
            yield open_file
    except (IOError, OSError) as e:
        logger.error(f"Error opening file {filename}: {e}")
        raise

def context_dnds(codon: str, mutated_codon: str) -> Dict[str, Any]:
    """
    Computes the dN/dS ratio for the given context.
    
    Args:
        codon: Original codon sequence
        mutated_codon: Mutated codon sequence
        
    Returns:
        Dictionary with N_sites, S_sites, and mutation_type
        
    Raises:
        ValueError: If codon is invalid
    """
    # Validate codon
    if len(codon) != 3 or not all(base in VALID_BASES for base in codon):
        raise ValueError(f"Invalid codon: {codon}. Codon must be a string of length 3 containing only A, C, G, T.")
    
    # Get the amino acid for the codon
    amino_acid = CODON_TO_AMINO.get(codon)
    if amino_acid is None:
        raise ValueError(f"Invalid codon: {codon}. Codon does not map to any amino acid.")
    
    # Count the number of synonymous and non-synonymous mutation sites
    N_sites = 0  # non-synonymous mutation sites
    S_sites = 0  # synonymous mutation sites
    
    for i in range(3):
        base = codon[i]
        for b in VALID_BASES:
            if b != base:
                new_codon = codon[:i] + b + codon[i + 1:]
                new_amino_acid = CODON_TO_AMINO.get(new_codon)
                if new_amino_acid is None:
                    continue
                
                if new_amino_acid == amino_acid:
                    S_sites += 1
                else:
                    N_sites += 1
    
    # Determine the type of mutation
    mutated_amino_acid = CODON_TO_AMINO.get(mutated_codon)
    if mutated_amino_acid == amino_acid:
        mutation_type = "synonymous"
    else:
        mutation_type = "non-synonymous"

    return {
        "N_sites": N_sites,
        "S_sites": S_sites,
        "mutation_type": mutation_type
    }

def context(    
    fasta: str, 
    vcf: str,
    model: str = "random", 
    alpha: float = 2.0, 
    beta: float = 1.0, 
    gamma: float = 1.0, 
    pi_a: float = 0.3,
    pi_c: float = 0.2,
    pi_g: float = 0.2,
    pi_t: float = 0.3,
    omega: float = 0.5) -> float:

    # Open the fasta file
    fasta_file = pysam.FastaFile(fasta)
    # Open the vcf file
    vcf_file = pysam.VariantFile(vcf)
    # Initialize the mutation model
    mutation_model = MutationModel(
        model_type=model,
        alpha=alpha,
        beta=beta,
        gamma=gamma,
        pi_a=pi_a,
        pi_c=pi_c,
        pi_g=pi_g,
        pi_t=pi_t,
        omega=omega
    )

    # Initialize the dN/dS ratio dict
    dnds_stats = {
        "N_sites": 0,  # non-synonymous mutation sites
        "S_sites": 0,  # synonymous mutation sites
        "synonymous": 0,
        "non_synonymous": 0
    }

    # Iterate by variant in the vcf file
    for record in vcf_file:
        chrom = record.chrom
        pos = record.pos
        ref_base = record.ref
        alt_base = record.alts[0] if record.alts else ref_base

        # random number 1-3
        seed = random.randint(1, 3)
        if seed == 1:
            codon = fasta_file.fetch(chrom, pos - 1, pos + 2).upper()
        elif seed == 2:
            codon = fasta_file.fetch(chrom, pos - 2, pos + 1).upper()
        else:
            codon = fasta_file.fetch(chrom, pos - 3, pos).upper()
        if len(codon) != 3 or not all(base in VALID_BASES for base in codon):
            logger.warning(f"Invalid codon at {chrom}:{pos} - skipping")
            continue
        # Get the amino acid for the codon
        amino_acid = CODON_TO_AMINO.get(codon)

        if amino_acid is None:
            logger.warning(f"Invalid codon at {chrom}:{pos} - skipping")
            continue
        # Calculate mutation probabilities for each base in the codon
        mutation = mutation_model.mutate(ref_base)
        if mutation is None:
            logger.warning(f"Mutation model returned None for base {ref_base} at {chrom}:{pos}")
            continue
        # Replace the reference base with the mutated base by the index of the seed
        mutated_codon = list(codon)
        mutated_codon[seed - 1] = mutation
        mutated_codon = "".join(mutated_codon)
        # Get the mutated amino acid
        mutated_amino_acid = CODON_TO_AMINO.get(mutated_codon)
        if mutated_amino_acid is None:
            logger.warning(f"Invalid mutated codon at {chrom}:{pos} - skipping")
            continue
        
        # Calculate dN/dS ratio
        dnds_codon_stat = context_dnds(codon, mutated_codon)
        dnds_stats["N_sites"] += dnds_codon_stat["N_sites"]
        dnds_stats["S_sites"] += dnds_codon_stat["S_sites"]
        if dnds_codon_stat["mutation_type"] == "synonymous":
            dnds_stats["synonymous"] += 1
        else:
            dnds_stats["non_synonymous"] += 1

    dnds = (dnds_stats["non_synonymous"] / dnds_stats["N_sites"]) / (dnds_stats["synonymous"] / dnds_stats["S_sites"])
    return dnds