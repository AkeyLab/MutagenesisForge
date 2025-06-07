from contextlib import contextmanager
import gzip
import pysam
import numpy as np
from collections import defaultdict

import yaml
import os

from .models import random_mutation, K2P, K3P, HKY85, JC69
from .utils import load_parameter_from_yaml, check_yaml_variable

"""
This module contains functions for creating a vcf file of random mutations.
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

def get_random_position_in_regions(regions: list):
    """
    Returns a random position within randomly selected region.

     Parameters:
        regions (list): list of regions in the bed file

    Returns:
        tuple: chromosome and position
    """
    region = regions[np.random.randint(0, len(regions))]
    region_start, region_end = map(int, region.split()[1:3])
    return region.split()[0], np.random.randint(region_start, region_end)


def get_base(fasta, chrom: str, pos: int):
    """
    Returns the reference base before and after provided chrom, position.
    """
    return fasta.fetch(chrom, pos - 1, pos).upper()


def is_random_pos_wanted(
    fasta: pysam.FastaFile, chrom: str, pos: str, before_base: str, after_base: str, ref_base: str, context_model: str
):
    """
    Determines whether the random position has same trinucleotide context.

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


def get_random_mut(before_base: str, after_base, ref_base, regions, fasta, context_model, model, alpha = None, beta = None, gamma = None, ratio = None):

    """
    Returns a random mutation in a random region that matches the specified criteria.

    Parameters:
        before_base (str): base before the position
        after_base (str): base after the position
        ref_base (str): reference base
        regions (list): list of regions in the bed file
        fasta (pysam.Fastafile): Fastafile object
        model (str): evolutionary model
        alpha (float): transition probability
        beta (float): transversion probability
        gamma (float): transversion probability

    Returns:
        tuple: chromosome, position, reference base, alternative base
    """

    # list of evolutionary models
    models = ["random", "K2P", "K3P"]
    # check if the model is valid
    if model not in models:
        raise ValueError(f"Model {model} is not valid. Choose from {models}")

    # initialize alpha, beta, and gamma for K2P and K3P models
    alpha, beta, gamma = None, None, None

    # add ratio for K2P and K3P models rather than using provided alpha and beta
    if model == "K2P":
        if ratio is None:
            raise ValueError("Transition-transversion ratio is required for K2P model")
        # mkae sure that ratio is in one of the following regex formats numeric/numeric or numeric:numeric
        if "/" in ratio:
            alpha, beta = map(float, ratio.split("/"))
        elif ":" in ratio:
            alpha, beta = map(float, ratio.split(":"))
        else:
            raise ValueError("Transition-transversion ratio is not in the correct format")
    
    if model == "K3P":
        if ratio is None:
            raise ValueError("Transition-transversion ratio is required for K3P model")
        # mkae sure that ratio is in one of the following regex formats numeric/numeric or numeric:numeric
        if "/" in ratio:
            alpha, beta, gamma = map(float, ratio.split("/"))
        elif ":" in ratio:
            alpha, beta, gamma = map(float, ratio.split(":"))
        else:
            raise ValueError("Transition-transversion ratio is not in the correct format")

    # get a random position from the regions
    is_wanted = False
    while not is_wanted:
        random_chr, random_pos = get_random_position_in_regions(regions)
        # determine whether the trinucleotide context is the same
        pos_base = get_base(fasta, random_chr, random_pos)
        # determine whether the trinucleotide context is the same
        is_wanted = is_random_pos_wanted(
            fasta, random_chr, random_pos, before_base, after_base, ref_base, context_model
        )

        # if the trinucleotide context is the same, get a random alternative allele
        if is_wanted:

            if model == "random":
                alt = random_mutation(ref_base)

            if model == "K2P":
                alt = K2P(ref_base, alpha, beta)
            
            if model == "K3P":
                alt = K3P(ref_base, alpha, beta, gamma)

            return random_chr, random_pos, ref_base, alt


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
