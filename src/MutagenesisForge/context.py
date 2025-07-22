from contextlib import contextmanager
import gzip
import pysam
from Bio import SeqIO
import numpy as np
from collections import defaultdict

import yaml
import os

from .mutation_model import MutationModel
from .utils import load_parameter_from_yaml, check_yaml_variable

"""
This module contains functions for creating a vcf file of random mutations.
"""

"""
Logic: Find random positions in the genome that match the created codon from a mutation following the mutation model.
Then, return the dN/dS ratio for the new positions
"""


@contextmanager
def my_open(filename: str, mode: str):
    """A wrapper for open/gzip.open logic as a context manager"""
    with (gzip.open(filename, mode + "t") if filename.endswith(".gz") else open(filename, mode)) as open_file:
        yield open_file


def get_trinucleotide_context(chrom: str, pos: str, fasta_file: pysam.Fastafile):
    """
    Gets the trinucleotide alleles at the specified position.

    Parameters:
        chrom (str): chromosome
        pos (int): position
        fasta_file (pysam.Fastafile): Fastafile object

    Returns:
        tuple: alleles at the position before, at, and after the specified position
    """

    return (
        fasta_file.fetch(chrom, pos - 2, pos - 1),
        fasta_file.fetch(chrom, pos - 1, pos),
        fasta_file.fetch(chrom, pos, pos + 1)
    )


def get_base(fasta, chrom: str, pos: int):
    """
    Returns the reference base before and after provided chrom, position.
    """
    return fasta.fetch(chrom, pos - 1, pos).upper()


def is_random_pos_wanted(
    fasta: pysam.FastaFile, 
    chrom: str, 
    pos: str, 
    before_base: str, 
    after_base: str, 
    ref_base: str, 
    context_model: str) -> bool:
    """
    Determines whether the random position has the correct trinucleotide context given the context model.

    Parameters:
        fasta (pysam.Fastafile): Fastafile object
        chrom (str): chromosome
        pos (int): position
        before_base (str): base before the position
        after_base (str): base after the position
        ref_base (str): reference base
        context_model (str): context model

    Returns:
        bool: True if the random position is wanted, False otherwise
    """
    pos_before = get_base(fasta, chrom, pos - 1)
    pos_after = get_base(fasta, chrom, pos + 1)
    pos_base = get_base(fasta, chrom, pos)
    # return True if the random position is wanted, False otherwise
    # wanted means the trinucleotide context is the same, and the before and after bases are the same
    
    # context model options

    # blind: no context model
    if context_model == "blind":
        return True
    # ra: reference allele context model
    if context_model == "ra":
        return (pos_base == ref_base)
    # ra_ba: reference allele and before allele context model
    if context_model == "ra_ba":
        return (pos_base == ref_base) and (pos_before == before_base)
    # ra_aa: reference allele and after allele context model
    if context_model == "ra_aa":
        return (pos_base == ref_base) and (pos_after == after_base)
    # codon: codon context model
    if context_model == "codon":
        return (pos_base == ref_base) and (pos_before == before_base) and (pos_after == after_base)
    else:
        raise ValueError(f"Context model {context_model} is not valid.")


# In-memory cache for context-matching positions
_context_match_cache = {}


def get_matching_positions(
    regions: list,
    fasta: pysam.FastaFile,
    before_base: str,
    after_base: str,
    ref_base: str,
    context_model: str,
) -> list[tuple[str, int]]:
    """
    Returns a list of all positions in regions that match the given context model.
    Uses in-memory cache to avoid redundant computation.
    """
    cache_key = (before_base, ref_base, after_base, context_model)

    if cache_key in _context_match_cache:
        return _context_match_cache[cache_key]

    matching = []
    for region in regions:
        chrom, start, end = region.split()[0], int(region.split()[1]), int(region.split()[2])
        for pos in range(start + 1, end - 1):
            try:
                if is_random_pos_wanted(fasta, chrom, pos, before_base, after_base, ref_base, context_model):
                    matching.append((chrom, pos))
            except Exception:
                continue  # skip over errors due to invalid sequence

    _context_match_cache[cache_key] = matching
    return matching


def get_random_mut(before_base: str, 
                   after_base: str, 
                   ref_base: str, 
                   regions: list, 
                   fasta: pysam.FastaFile, 
                   context_model: str, 
                   model: MutationModel, 
                   ) -> tuple:

    """
    Returns a random mutation in a random region that matches the specified criteria.

    Parameters:
        before_base (str): base before the position
        after_base (str): base after the position
        ref_base (str): reference base
        regions (list): list of regions in the bed file
        fasta (pysam.Fastafile): Fastafile object
        context_model (str): context model
        model (str): mutation model type
        alpha (float): alpha parameter for mutation model
        beta (float): beta parameter for mutation model
        gamma (float): gamma parameter for mutation model
        pi_a (float): frequency of A base in mutation model
        pi_c (float): frequency of C base in mutation model
        pi_g (float): frequency of G base in mutation model
        pi_t (float): frequency of T base in mutation model

    Returns:
        tuple: chromosome, position, reference base, alternative base
    """
    
    matching_positions = get_matching_positions(
        regions,
        fasta,
        before_base,
        after_base,
        ref_base,
        context_model
    )

    if not matching_positions:
        raise ValueError(f"No matching context found for {before_base}-{ref_base}-{after_base} under model {context_model}")

    # randomly select a position from the matching positions
    random_chr, random_pos = matching_positions[np.random.randint(len(matching_positions))]
    
    # Mutate the base
    alt = model.mutate(ref_base)

    
    return random_chr, random_pos, ref_base, alt


def context_dnds(codon, mutated_codon) -> dict:
    """
    Computes the dN/dS ratio for the given context.
    """
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

    # check if the codon is valid
    if len(codon) != 3 or not all(base in bases for base in codon):
        raise ValueError(f"Invalid codon: {codon}. Codon must be a string of length 3 containing only A, C, G, T.")
    # get the amino acid for the codon
    amino_acid = codon_to_amino.get(codon, None)
    if amino_acid is None:
        raise ValueError(f"Invalid codon: {codon}. Codon does not map to any amino acid.")
    # count the number of synonymous and non-synonymou mutation cites
    N_sites = 0  # non-synonymous mutation sites
    S_sites = 0  # synonymous mutations sites
    for i in range(3):
        # get the base at the current position
        base = codon[i]
        # iterate through the bases
        for b in bases:
            if b != base:
                # create a new codon with the mutated base
                new_codon = codon[:i] + b + codon[i + 1:]
                # get the amino acid for the new codon
                new_amino_acid = codon_to_amino.get(new_codon, None)
                if new_amino_acid is None:
                    continue
                # check if the new amino acid is the same as the original amino acid
                if new_amino_acid == amino_acid:
                    S_sites += 1
                else:
                    N_sites += 1
    # determine the type of mutation
    if codon_to_amino.get(mutated_codon, None) == amino_acid:
        mutation_type = "synonymous"
    else:
        mutation_type = "non-synonymous"

    return {
        "N_sites": N_sites,
        "S_sites": S_sites,
        "mutation_type": mutation_type
    }


def context(fasta: str, 
            vcf: str, 
            bed: str, 
            model:str = "random", 
            alpha:float = None, 
            beta:float = None, 
            gamma:float = None, 
            pi_a:float = None,
            pi_c:float = None,
            pi_g:float = None,
            pi_t:float = None,
            context_model:str = "codon"):
    """
    Run the simulation and return calculate the dN/dS ratio.
    """
    # read in the fasta file
    fasta = pysam.Fastafile(fasta)
    # read in the vcf file
    vcf_file = pysam.VariantFile(vcf)
    # read in the bed file
    regions = []

    if bed is not None:
        with my_open(bed, "r") as f:
            for line in f:
                line = line.strip()
                if not line.startswith("#"):
                    regions.append(line)
    else:
        with my_open(fasta, "r") as f:
            for record in SeqIO.parse(f, "fasta"):
                chrom = record.id
                length = len(record.seq)
                regions.append(f"{chrom}\t0\t{length}")
    
    # create a mutation model
    mutation_model = MutationModel(
        model_type=model,
        gamma=gamma,
        alpha=alpha,
        beta=beta,
        pi_a=pi_a,
        pi_c=pi_c,
        pi_g=pi_g,
        pi_t=pi_t
    )
    # create a dictionary to store the dN/dS data
    dnds_data = defaultdict(lambda: {"N_sites": 0, "S_sites": 0, "synonymous": 0, "non_synonymous": 0})
    # iterate through the vcf file
    for record in vcf_file:
        # get the chromosome and position
        chrom = record.chrom
        pos = record.pos
        # get the reference base and alternative base
        ref_base = record.ref
        alt_base = record.alts[0] if record.alts else None
        if alt_base is None:
            alt_base = ref_base  # if no alt base, use ref base
            continue
        # get the trinucleotide context
        before_base, ref_base, after_base = get_trinucleotide_context(chrom, pos, fasta)

        # find matching random matching position within the regions
        matching_positions = get_matching_positions(
            regions,
            fasta,
            before_base,
            after_base,
            ref_base,
            context_model
        )
        # select a random position from the matching positions
        random_chr, random_pos = matching_positions[np.random.randint(len(matching_positions))]
        # get codon context
        codon = get_trinucleotide_context(random_chr, random_pos, fasta)
        # get the mutated codon
        random_alt_base = mutation_model.mutate(ref_base)
        mutated_codon = codon[0] + random_alt_base + codon[2]  # replace the middle base with the mutated base

        # calculate dN/dS statistics for the mutated codon
        dnds_stats = context_dnds(codon, mutated_codon, model)
        # update the dN/dS data
        dnds_data["N_sites"] += dnds_stats["N_sites"]
        dnds_data["S_sites"] += dnds_stats["S_sites"]
        if dnds_stats["mutation_type"] == "synonymous":
            dnds_data["synonymous"] += 1
        else:
            dnds_data["non_synonymous"] += 1
    
    # calculate the dN/dS ratio
    dN = dnds_data["non_synonymous"] / dnds_data["N_sites"] if dnds_data["N_sites"] > 0 else 0
    dS = dnds_data["synonymous"] / dnds_data["S_sites"] if dnds_data["S_sites"] > 0 else 0
    dnds_ratio = dN / dS if dS > 0 else float('inf')  # handle division by zero
    # return the dN/dS ratio and the dN and dS values
    """
    return {
        "dN": dN,
        "dS": dS,
        "dN/dS": dnds_ratio,
        "N_sites": dnds_data["N_sites"],
        "S_sites": dnds_data["S_sites"],
        "synonymous": dnds_data["synonymous"],
        "non_synonymous": dnds_data["non_synonymous"]
    }
    """
    return dnds_ratio

def create_vcf_file(input_file, output_file):
    """
    Create a vcf file from the input file.

    Parameters:
        input_file (str): input file path
        output_file (str): output vcf file path

    Returns:
        None (writes to output vcf file)
    """
    # read in the variants from the input file
    with open(input_file, "r") as f:
        variants = f.readlines()
    # create the vcf header
    vcf_header = "##fileformat=VCFv4.3\n"
    vcf_header += '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    vcf_content = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
    # create a dictionary to store the variants by chromosome and position
    variant_dict = defaultdict(dict)
    # iterate through the variants and create the vcf content
    for variant in variants:
        variant_data = variant.strip().split("\t")
        if variant_data[0] == "CHR":
            continue
        chrom = int(variant_data[0])
        pos = int(variant_data[1])
        ref = variant_data[2]
        alt = variant_data[-1]
        vcf_line = f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\tGT\t0/1\n"
        variant_dict[chrom][pos] = vcf_line
        vcf_content += vcf_line
    # write the vcf content to the output file
    with open(output_file, "w") as f:
        f.write(vcf_header)
        # write the vcf content: sort the variants by chromosome and position
        for chrom in sorted(variant_dict.keys()):
            for pos in sorted(variant_dict[chrom].keys()):
                f.write(variant_dict[chrom][pos])


def indy_vep(vep_string, num, output):
    """
    Modift vep call output file to individualize.
    """
    return vep_string.replace("-o", f" -o {output}{num}.vep ")


def vcf_constr(bed_file, mut_file, fasta_file, output,
                sim_num, vep_call,
                model = "random", alpha = None, beta = None, gamma = None, context_model = "codon"):
    """
    Create a vcf file of random mutations given a bed file, mutation file, fasta file, output
    file, transition-transversion ratio, number of simulations, and whether to run vep call.

    Parameters:
        bed_file (str): path to bed file
        mut_file (str): path to mutation file
        fasta_file (str): path to fasta file
        output (str): output file name
        tstv (float): transition-transversion ratio
        sim_num (int): number of simulations
        vep_call (bool): run vep call

    Returns:
        None (creates a vcf file of random mutations)
    """
    # convert fasta file path into fastafile object
    print(f"cotext model: {context_model}")
    print(f"model: {model}")
    fasta = pysam.Fastafile(fasta_file)

    # read in regions from bed file
    regions = []
    with my_open(bed_file, "r") as f:
        for line in f:
            line = line.strip()
            regions.append(line)

    for i in range(sim_num):
        # code goes here
        # create output file names based on the iteration
        file_shell = output + str(i) + ".txt"
        vcf_shell = output + str(i) + ".vcf"
        with my_open(mut_file, "r") as f, my_open(file_shell, "w") as o:
            header = f.readline().strip().split()
            header_dict = dict(zip(header, range(len(header))))
            chr_pos_dict = {}
            # count = 0
            for line in f:
                if line.startswith("#"):
                    continue
                line = line.strip().split()
                chromosome = line[0]
                position = line[1]
                before_base, ref_base, after_base = get_trinucleotide_context(
                    str(chromosome), int(position), fasta
                )

                # find randomized mutations
                add_one_random_mut = False
                while not add_one_random_mut:

                    random_chr, random_pos, ref_base, alt = get_random_mut(
                        before_base, after_base, ref_base, regions, fasta, context_model, model, alpha, beta, gamma
                    )
                    chr_pos = random_chr + "_" + str(random_pos)
                    if chr_pos not in chr_pos_dict:
                        chr_pos_dict[chr_pos] = 1
                        add_one_random_mut = True
                        out_line = [
                            random_chr,
                            str(random_pos),
                            ref_base,
                            before_base,
                            after_base,
                            alt,
                        ]
                        o.write("\t".join([str(x) for x in out_line]) + "\n")
        create_vcf_file(file_shell, vcf_shell)

        # flag for vep run
        if vep_call:

            if not check_yaml_variable("parameters.yaml", "vep_tool_path"):
                print("vep_tool_path not found in parameters.yaml")
                return
            # run vep call on created vcf file
            # current path is just the path for my personal computer
            """
            vep = indy_vep(
                load_parameter_from_yaml("parameters.yaml", "vep_tool_path"),
                i,
                output,
            )
            """
            vep = load_parameter_from_yaml("parameters.yaml", "vep_tool_path") + " -i " + vcf_shell + " -o " + output + str(sim_num)
            print(vep)
            os.system(vep)


if __name__ == "__main__":
    pass
