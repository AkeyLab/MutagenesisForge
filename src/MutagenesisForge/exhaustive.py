import collections
import pysam
import numpy as np
from typing import Optional

from .mutation_model import MutationModel 

"""
This module contains evolutionary models for simulating mutations
Now supports BED file input for extracting coding regions.
"""

def exhaustive(fasta: str, 
               bed: Optional[str] = None, 
               by_read: bool = False, 
               model:str = "random", 
               alpha:str = None, 
               beta:str = None, 
               gamma:str = None,
               pi_a:str = None,
               pi_c:str = None,
               pi_g:str = None,
               pi_t:str = None):
    
    codon_to_amino = {
        "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
        "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M", "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
        "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
        "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T", "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
        "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*", "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
        "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K", "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
        "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W", "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
        "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R", "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G"
    }

    bases = {"A", "C", "G", "T"}

    # intialize mutation model
    if model == "random":
        mutation_model = MutationModel(
            model_type="random"
            )
    elif model == "JC69":
        if gamma is None:
            raise ValueError("Gamma parameter must be provided for JC69 model.") 
        mutation_model = MutationModel(
            model_type="JC69", 
            gamma=float(gamma))
    elif model == "K2P":
        if alpha is None or beta is None:
            raise ValueError("Alpha and beta parameters must be provided for K2P model.")
        mutation_model = MutationModel(model_type="K2P", 
                                       gamma=float(gamma), 
                                       alpha=float(alpha), 
                                       beta=float(beta))
    elif model == "F81":
        if pi_a is None or pi_c is None or pi_g is None or pi_t is None:
            raise ValueError("All base frequencies (pi_a, pi_c, pi_g, pi_t) must be provided for F81 model.")
        mutation_model = MutationModel(model_type="F81", 
                                       gamma=float(gamma), 
                                       pi_a=float(pi_a), 
                                       pi_c=float(pi_c), 
                                       pi_g=float(pi_g), 
                                       pi_t=float(pi_t))
    elif model == "HKY85":
        if alpha is None or beta is None or pi_a is None or pi_c is None or pi_g is None or pi_t is None:
            raise ValueError("Alpha, beta, and all base frequencies (pi_a, pi_c, pi_g, pi_t) must be provided for HKY85 model.")
        mutation_model = MutationModel(model_type="HKY85", 
                                       gamma=float(gamma), 
                                       alpha=float(alpha), 
                                       beta=float(beta), 
                                       pi_a=float(pi_a), 
                                       pi_c=float(pi_c), 
                                       pi_g=float(pi_g), 
                                       pi_t=float(pi_t))

    synonymous_muts_per_codon = collections.defaultdict(int)
    missense_muts_per_codon = collections.defaultdict(int)
    nonsense_muts_per_codon = collections.defaultdict(int)

    for codon in codon_to_amino:
        for pos in [0, 1, 2]:
            base = codon[pos]
            amino = codon_to_amino[codon]
            possible_mutations = bases - {base}

            for mutated_base in possible_mutations:
                # Calculate mutation probability using the mutation model
                prob = mutation_model.get_mutation_probability(base, mutated_base)
                mutated_codon = codon[:pos] + mutated_base + codon[pos + 1:]
                mutated_amino = codon_to_amino.get(mutated_codon, None)

                if mutated_amino is None:
                    continue

                if mutated_amino == "*":
                    nonsense_muts_per_codon[codon] += prob
                elif amino != mutated_amino:
                    missense_muts_per_codon[codon] += prob
                else:
                    synonymous_muts_per_codon[codon] += prob

    fasta = pysam.FastaFile(fasta)

    data = {
        "has_ATG_start": [],
        "stop_codon_end": [],
        "cds_length": [],
        "num_codons": [],
        "num_synonymous": [],
        "num_missense": [],
        "num_nonsense": [],
    }

    if bed:
        with open(bed) as bed:
            regions = [line.strip().split() for line in bed if not line.startswith("#")]
    else:
        regions = [(ref, 0, fasta.get_reference_length(ref)) for ref in fasta.references]

    for ref, start, end in regions:
        start, end = int(start), int(end)
        seq = fasta.fetch(ref, start, end)

        has_atg_start = seq[:3] == "ATG"
        stop_codons = ["TAA", "TGA", "TAG"]
        stop_codon_end = seq[-3:] in stop_codons

        codon_frequency = collections.Counter(seq[i:i + 3] for i in range(0, len(seq), 3))
        codon_frequency = {
            codon: count for codon, count in codon_frequency.items() if len(codon) == 3
        }

        num_synonymous = sum(
            synonymous_muts_per_codon.get(codon, 0) * count
            for codon, count in codon_frequency.items()
        )
        num_missense = sum(
            missense_muts_per_codon.get(codon, 0) * count
            for codon, count in codon_frequency.items()
        )
        num_nonsense = sum(
            nonsense_muts_per_codon.get(codon, 0) * count
            for codon, count in codon_frequency.items()
        )

        data["has_ATG_start"].append(has_atg_start)
        data["stop_codon_end"].append(stop_codon_end)
        data["num_codons"].append(sum(codon_frequency.values()))
        data["num_synonymous"].append(num_synonymous)
        data["num_missense"].append(num_missense)
        data["num_nonsense"].append(num_nonsense)

    fasta.close()

    dnds_method_1 = (sum(data["num_missense"]) + sum(data["num_nonsense"])) / max(sum(data["num_synonymous"]), 1)
    dnds_method_2 = []
    for i in range(len(data["has_ATG_start"])):
        if data["num_synonymous"][i] == 0:
            continue
        dnds_method_2.append(
            (data["num_missense"][i] + data["num_nonsense"][i]) / data["num_synonymous"][i]
        )

    dnds2_mean = np.mean(dnds_method_2) if dnds_method_2 else float('nan')

    return dnds2_mean if by_read else dnds_method_1