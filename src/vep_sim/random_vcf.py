from contextlib import contextmanager
import gzip
import pysam
import numpy as np
from collections import defaultdict

import yaml
import os

from .utils import load_parameter_from_yaml, check_yaml_variable

"""
This module contains functions for creating a vcf file of random mutations.
"""


@contextmanager
def my_open(filename: str, mode: str):
    """A wrapper for open/gzip.open logic as a context manager"""
    if filename.endswith(".gz"):
        open_file = gzip.open(filename, mode + "t")
    else:
        open_file = open(filename, mode)
    yield open_file
    open_file.close()


def get_trinucleotide_context(chrom, pos, fasta_file):
    """
    Gets the trinucleotide alleles at the specified position.

    Parameters:
        chrom (str): chromosome
        pos (int): position
        fasta_file (pysam.Fastafile): Fastafile object

    Returns:
        tuple: alleles at the position before, at, and after the specified position
    """

    # Fetch allele at the specified position
    allele_at_position = fasta_file.fetch(chrom, pos - 1, pos)

    # Fetch allele at the position before
    allele_before_position = fasta_file.fetch(chrom, pos - 2, pos - 1)

    # Fetch allele at the position after
    allele_after_position = fasta_file.fetch(chrom, pos, pos + 1)

    return allele_before_position, allele_at_position, allele_after_position


def get_base(fasta, chrom: str, position: str):
    """
    Returns the reference base before and after provided chrom, position.
    """
    pos = int(position)
    base = fasta.fetch(chrom, pos - 1, pos).upper()
    return base


def get_random_position_in_regions(regions):
    """
    Returns a random position in the regions provided in bed file.

    Parameters:
        regions (list): list of regions in the bed file

    Returns:
        tuple: chromosome and position
    """
    # region is a string, e.g. "chr1 100 200"
    # return a random position in the region
    region = regions[np.random.randint(0, len(regions))]
    region_start = int(region.split()[1])
    region_end = int(region.split()[2])
    pos = np.random.randint(region_start, region_end)
    return region.split()[0], pos


def is_random_pos_wanted(
    fasta, random_chr, random_pos, before_base, after_base, ref_base
):
    """
    Determines whether the random position has same trinucleotide context.

    Parameters:
        fasta (pysam.Fastafile): Fastafile object
        random_chr (str): chromosome
        random_pos (int): position
        before_base (str): base before the position
        after_base (str): base after the position
        ref_base (str): reference base

    Returns:
        bool: True if the random position is wanted, False otherwise
    """
    # return True if the random position is wanted, False otherwise
    # wanted means the trinucleotide context is the same, and the before and after bases are the same
    pos_base = get_base(fasta, random_chr, random_pos)
    if pos_base != ref_base:
        return False
    pos_before_base = get_base(fasta, random_chr, random_pos - 1)
    if pos_before_base != before_base:
        return False
    pos_after_base = get_base(fasta, random_chr, random_pos + 1)
    if pos_after_base != after_base:
        return False
    return True


def get_random_mut(before_base, after_base, ref_base, regions, fasta, tstv):
    """
    Returns a random mutation in a random region that matches the specified criteria.

    Parameters:
        before_base (str): base before the position
        after_base (str): base after the position
        ref_base (str): reference base
        regions (list): list of regions in the bed file
        fasta (pysam.Fastafile): Fastafile object
        tstv (float): transition-transversion ratio

    Returns:
        tuple: chromosome, position, reference base, alternative base
    """
    # get a random position from the regions
    is_wanted = False
    while not is_wanted:
        random_chr, random_pos = get_random_position_in_regions(regions)
        # determine whether the trinucleotide context is the same
        pos_base = get_base(fasta, random_chr, random_pos)
        # determine whether the trinucleotide context is the same
        is_wanted = is_random_pos_wanted(
            fasta, random_chr, random_pos, before_base, after_base, ref_base
        )
        # if the trinucleotide context is the same, get a random alternative allele
        if is_wanted:
            # p(s) = tstv/(tstv+2)
            # implementation of tstv ratio
            if pos_base == "A":
                alt = str(
                    np.random.choice(
                        ["G", "T", "C"],
                        p=[tstv / (tstv + 2), 1 / (tstv + 2), 1 / (tstv + 2)],
                    )
                )
            if pos_base == "G":
                alt = str(
                    np.random.choice(
                        ["A", "T", "C"],
                        p=[tstv / (tstv + 2), 1 / (tstv + 2), 1 / (tstv + 2)],
                    )
                )
            if pos_base == "T":
                alt = str(
                    np.random.choice(
                        ["C", "A", "G"],
                        p=[tstv / (tstv + 2), 1 / (tstv + 2), 1 / (tstv + 2)],
                    )
                )
            if pos_base == "C":
                alt = str(
                    np.random.choice(
                        ["T", "A", "G"],
                        p=[tstv / (tstv + 2), 1 / (tstv + 2), 1 / (tstv + 2)],
                    )
                )
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
        #print(vcf_line) debugging remnant
        variant_dict[chrom][pos] = vcf_line
        vcf_content += vcf_line
    # write the vcf content to the output file
    with open(output_file, "w") as f:
        f.write(vcf_header)
        # write the vcf content: sort the variants by chromosome and position
        for chrom in sorted(variant_dict.keys()):
            for pos in sorted(variant_dict[chrom].keys()):
                f.write(variant_dict[chrom][pos])


# Removed the following code from the function vcf_constr:
#     yaml finding because it is in utils


def indy_vep(vep_string, num, output):
    """
    Individualizes each vep call output file.
    """
    motif = "-o"
    motif_index = vep_string.find(motif)
    return (
        vep_string[:motif_index]
        + " "
        + output
        + str(num)
        + " .vep"
        + vep_string[motif_index:]
    )


def vcf_constr(bed_file, mut_file, fasta_file, output, tstv, sim_num, vep_call):
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
                        before_base, after_base, ref_base, regions, fasta, tstv
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

            if not check_yaml_variable("../../parameters.yaml", "vep_tool_path"):
                print("vep_tool_path not found in parameters.yaml")
                return
            # run vep call on created vcf file
            # current path is just the path for my personal computer
            vep = indy_vep(
                load_parameter_from_yaml("../../parameters.yaml", "vep_tool_path"),
                i,
                output,
            )
            os.system(vep)


if __name__ == "__main__":
    pass
