from contextlib import contextmanager
import gzip
import numpy as np

fasta = "test.fa"

# think about implementing function such as get_mutation for Mutation class addition
# this will help with having each muation class act as a line from within the vcf


class Mutation:
    def __init__(self, line):
        self.chrom = line.split()[0]
        self.pos = line.split()[1]
        self.ref = line.split()[3]
        self.alt = line.split()[4]

    def tri_nuc_context(self):
        # Get the allele at the specified position
        allele_at_position = fasta.fetch(self.chrom, self.pos - 1, self.pos)

        # Get the allele at the position before
        allele_before_position = fasta.fetch(self.chrom, self.pos - 2, self.pos - 1)

        # Get the allele at the position after
        allele_after_position = fasta.fetch(self.chrom, self.pos, self.pos + 1)

        return [allele_before_position, allele_at_position, allele_after_position]
    
    def new_mut(self):
        self.tri_nuc_context()

class VCF:
    def __init__(self, filename):
        self.filename = filename
        self.mutations = []
        self.parse_vcf()


