from contextlib import contextmanager
import gzip
import pysam
import numpy as np
from collections import defaultdict

import yaml
import os

# return vcf file of random mutations

# for vep call, use a combination of yaml input and user input for the files run through vep
# ex. yaml input then remove -o and -i and replace with user input through click via main.py
# need to include more input flags potentially

#working
@contextmanager
def my_open(filename: str, mode: str):
    '''A wrapper for open/gzip.open logic as a context manager'''
    if filename.endswith('.gz'):
        open_file = gzip.open(filename, mode+'t')
    else:
        open_file = open(filename, mode)
    yield open_file
    open_file.close()

# my_open test
#with my_open("chromosome10_test.vcf", "r") as f:
#    for line in f:
#        if line.startswith('#'):
#            continue
            #print(line)

# chromosome must be string
#is working
def get_trinucleotide_context(chrom, pos, fasta_file):
    # Get the allele at the specified position
    allele_at_position = fasta_file.fetch(chrom, pos - 1, pos)

    # Get the allele at the position before
    allele_before_position = fasta_file.fetch(chrom, pos - 2, pos - 1)

    # Get the allele at the position after
    allele_after_position = fasta_file.fetch(chrom, pos, pos + 1)

    return allele_before_position, allele_at_position, allele_after_position

fasta = '/projects/AKEY/akey_vol2/References/Genomes/hs37d5/hs37d5.fa'
fasta_file = pysam.Fastafile(fasta)
print(get_trinucleotide_context('12', 109845669, fasta_file))


#get_trinucleotide_context("chr1", 100, pysam.Fastafile("test.fa"))

#working
def get_base(fasta, chrom: str, position: str):
    '''
    Return reference base before and after provided chrom, position
    '''
    pos = int(position)
    # want surrounding two bases, fasta is zero indexed while position is unit
    base = fasta.fetch(chrom, pos-1, pos).upper()
    return base
print(get_base(fasta_file, '12', '109845669'))

#is working
def get_random_position_in_regions(regions):
    #region is a string, e.g. "chr1 100 200"
    #return a random position in the region
    region = regions[np.random.randint(0, len(regions))]
    region_start = int(region.split()[1])
    region_end = int(region.split()[2])
    pos = np.random.randint(region_start, region_end)
    return region.split()[0], pos

input_bed_file = '/projects/AKEY/akey_vol2/huixinx/Projects/01.eGTEx/NWGC/04.fig3/02.exp_mis_to_syn_ratio/step12.problematic.bed'

regions = []
with my_open(input_bed_file, 'r') as f:
    for line in f:    
        line = line.strip()
        regions.append(line)
print(get_random_position_in_regions(regions))


def is_random_pos_wanted(fasta, random_chr, random_pos, before_base, after_base, ref_base):
    #return True if the random position is wanted, False otherwise
    #wanted means the trinucleotide context is the same, and the before and after bases are the same
    pos_base = get_base(fasta, random_chr, random_pos)
    if pos_base != ref_base:
        return False
    pos_before_base = get_base(fasta, random_chr, random_pos-1)
    if pos_before_base != before_base:
        return False
    pos_after_base = get_base(fasta, random_chr, random_pos+1)
    if pos_after_base != after_base:
        return False
    return True


def get_random_mut(before_base, after_base, ref_base, regions, fasta, tstv):
    #before_base and after_base are strings of length 1
    #ref is a string of length 3
    #regions is an array of strings, each string is a region in the bed file
    #return a random mutation in a random region, and the mutation should match the following criteria: same trinucleotide context, random position from the bed regions, random alternative allele
    #return a tuple of (chr, pos, ref, alt)
    #proabilities account for transverstion-transition bias

    # get a random position from the regions
    is_wanted = False
    while not is_wanted:
        random_chr, random_pos = get_random_position_in_regions(regions)
        # determine whether the trinucleotide context is the same
        pos_base = get_base(fasta, random_chr, random_pos)
        # determine whether the trinucleotide context is the same
        is_wanted = is_random_pos_wanted(fasta, random_chr, random_pos, before_base, after_base, ref_base)
        # if the trinucleotide context is the same, get a random alternative allele
        if is_wanted:
            #tstv ratio will be added
            #p(s) = tstv/(tstv+2)
            if pos_base == 'A':
                alt = str(np.random.choice(['G', 'T', 'C'], p=[tstv/(tstv+2), 1/(tstv+2), 1/(tstv+2)]))
            if pos_base == 'G':
                alt = str(np.random.choice(['A', 'T', 'C'], p=[tstv/(tstv+2), 1/(tstv+2), 1/(tstv+2)]))
            if pos_base == 'T':
                alt = str(np.random.choice(['C', 'A', 'G'], p=[tstv/(tstv+2), 1/(tstv+2), 1/(tstv+2)]))
            if pos_base == 'C':
                alt = str(np.random.choice(['T', 'A', 'G'], p=[tstv/(tstv+2), 1/(tstv+2), 1/(tstv+2)]))
            return random_chr, random_pos, ref_base, alt

print(get_random_mut('A', 'C', 'T', regions, fasta_file, 2.0))


def create_vcf_file(input_file, output_file):

    with open(input_file, 'r') as f:
        variants = f.readlines()

    vcf_header = '##fileformat=VCFv4.3\n'
    vcf_header += '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    vcf_content = '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n'

    variant_dict = defaultdict(dict)

    for variant in variants:
        variant_data = variant.strip().split('\t')
        if variant_data[0]=="CHR":
            continue
        chrom = int(variant_data[0])
        pos = int(variant_data[1])
        ref = variant_data[2]
        # the last element of variant_data is the alternative allele
        alt = variant_data[-1]

        vcf_line = f'{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\tGT\t0/1\n'
        print(vcf_line)
        variant_dict[chrom][pos] = vcf_line
        vcf_content += vcf_line

    with open(output_file, 'w') as f:
        f.write(vcf_header)
        # write the vcf content: sort the variants by chromosome and position
        for chrom in sorted(variant_dict.keys()):
            for pos in sorted(variant_dict[chrom].keys()):
                f.write(variant_dict[chrom][pos])


def load_parameter_from_yaml(file_path, parameter_name):
    with open(file_path, 'r') as file:
        params = yaml.safe_load(file)
        if parameter_name in params:
            return params[parameter_name]
        else:
            print(f"Parameter '{parameter_name}' not found in {file_path}")
            return None


def vcf_constr(bed_file, mut_file, fasta_file, output, tstv, sim_num, vep_call = False):
    # convert fasta file path into fastafile object
    fasta = pysam.Fastafile(fasta_file)
    
    # read in regions from bed file
    regions = []
    with my_open(bed_file, 'r') as f:
        for line in f:
            line = line.strip()
            regions.append(line)

    for i in range(sim_num):
        # create output file names based on the iteration
        file_shell = output + str(i) + '.txt'
        vcf_shell = output + str(i) + '.vcf'
        with my_open(mut_file, 'r') as f, my_open(file_shell, 'w') as o:
            header = f.readline().strip().split()
            header_dict = dict(zip(header, range(len(header))))
            chr_pos_dict = {}
            # count = 0
            for line in f:
                ''' this is for testing purposes
                count += 1
                if count > 30:
                    break
                    '''
                if line.startswith('#'):
                    continue
                line = line.strip().split()
                chromosome = line[0]
                position = line[1]
                before_base, ref_base, after_base = get_trinucleotide_context(str(chromosome), int(position), fasta)

                # find randomized mutations
                add_one_random_mut = False
                while not add_one_random_mut:
                    random_chr, random_pos, ref_base, alt = get_random_mut(before_base, after_base, ref_base, regions, fasta, tstv)
                    chr_pos = random_chr + "_" + str(random_pos)
                    if chr_pos not in chr_pos_dict:
                        chr_pos_dict[chr_pos] = 1
                        add_one_random_mut = True
                        out_line = [random_chr, str(random_pos), ref_base, before_base, after_base, alt]
                        o.write('\t'.join([str(x) for x in out_line]) + '\n')
        create_vcf_file(file_shell, vcf_shell)
        
        # flag for vep run
        if vep_call:
            # run vep call on created vcf file
            # current path is just the path for my personal computer
            vep = load_parameter_from_yaml('/projects/AKEY/akey_vol2/cooper/vep-sim/parameters.yaml', 'vep_tool_path')
            os.system(vep_call)
    

vcf_constr('/projects/AKEY/akey_vol2/huixinx/Projects/01.eGTEx/NWGC/04.fig3/02.exp_mis_to_syn_ratio/step12.problematic.bed', 
        '/projects/AKEY/akey_vol2/huixinx/Projects/01.eGTEx/NWGC/04.fig3/02.exp_mis_to_syn_ratio/docker_stringent.nwgc.rep2.raw_bb_p_lt_10_8.filtered10.txt',
        '/projects/AKEY/akey_vol2/References/Genomes/hs37d5/hs37d5.fa',
        'output.out')