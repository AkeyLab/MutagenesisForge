# using UCSC or 1000 genome database (hg38) create pseudo-codons for simulated mutations

from pyfaidx import Fasta

# Path to your hg38 reference genome FASTA file
fasta_file = 'path_to_your_hg38_reference.fasta'

def get_bases_around_location(chromosome, position):
    # Load the reference genome FASTA file
    genome = Fasta(fasta_file)

    # Retrieve the sequence from the reference genome
    base_before = genome[chromosome][position - 1].seq
    base = genome[chromosome][position].seq
    base_after = genome[chromosome][position + 1].seq
    
    return base_before, base, base_after


