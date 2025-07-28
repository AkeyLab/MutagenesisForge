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
DEFAULT_MAX_MATCHES = 1000

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


def get_trinucleotide_context(chrom: str, pos: int, fasta_file: pysam.Fastafile) -> Optional[Tuple[str, str, str]]:
    """
    Gets the trinucleotide alleles at the specified position.

    Args:
        chrom: chromosome
        pos: position (1-based)
        fasta_file: Fastafile object

    Returns:
        Tuple of (before_base, current_base, after_base) or None if invalid bases found
    """
    try:
        # Check bounds - need at least position 2 and not at chromosome end
        chrom_length = fasta_file.get_reference_length(chrom)
        if pos < 2 or pos >= chrom_length:
            logger.debug(f"Position {chrom}:{pos} too close to chromosome boundary (length: {chrom_length})")
            return None
            
        before = fasta_file.fetch(chrom, pos - 2, pos - 1).upper()
        current = fasta_file.fetch(chrom, pos - 1, pos).upper()
        after = fasta_file.fetch(chrom, pos, pos + 1).upper()
        
        # Skip positions with ambiguous bases or empty sequences
        if not all(base and base in VALID_BASES for base in [before, current, after]):
            logger.debug(f"Skipping position {chrom}:{pos} due to ambiguous bases: {before}-{current}-{after}")
            return None
            
        return (before, current, after)
    except Exception as e:
        logger.warning(f"Error getting trinucleotide context at {chrom}:{pos}: {e}")
        return None


def get_base(fasta: pysam.FastaFile, chrom: str, pos: int) -> Optional[str]:
    """
    Returns the reference base at the provided chromosome and position.
    
    Args:
        fasta: FastaFile object
        chrom: chromosome
        pos: position (1-based)
        
    Returns:
        Base at position or None if invalid
    """
    try:
        # Check bounds
        chrom_length = fasta.get_reference_length(chrom)
        if pos < 1 or pos > chrom_length:
            logger.debug(f"Position {chrom}:{pos} out of bounds (length: {chrom_length})")
            return None
            
        base = fasta.fetch(chrom, pos - 1, pos).upper()
        return base if base and base in VALID_BASES else None
    except Exception as e:
        logger.warning(f"Error getting base at {chrom}:{pos}: {e}")
        return None


def is_random_pos_wanted(
    fasta: pysam.FastaFile, 
    chrom: str, 
    pos: int, 
    before_base: str, 
    after_base: str, 
    ref_base: str, 
    context_model: str
) -> bool:
    """
    Determines whether the random position has the correct trinucleotide context given the context model.

    Args:
        fasta: Fastafile object
        chrom: chromosome
        pos: position
        before_base: base before the position
        after_base: base after the position
        ref_base: reference base
        context_model: context model

    Returns:
        True if the random position is wanted, False otherwise
    """
    pos_before = get_base(fasta, chrom, pos - 1)
    pos_after = get_base(fasta, chrom, pos + 1)
    pos_base = get_base(fasta, chrom, pos)
    
    # Skip if any base is invalid
    if any(base is None for base in [pos_before, pos_after, pos_base]):
        return False
    
    # Context model options
    if context_model == "blind":
        return True
    elif context_model == "ra":
        return pos_base == ref_base
    elif context_model == "ra_ba":
        return (pos_base == ref_base) and (pos_before == before_base)
    elif context_model == "ra_aa":
        return (pos_base == ref_base) and (pos_after == after_base)
    elif context_model == "codon":
        return (pos_base == ref_base) and (pos_before == before_base) and (pos_after == after_base)
    else:
        raise ValueError(f"Context model {context_model} is not valid.")


@lru_cache(maxsize=500)
def get_matching_positions_cached(
    regions_tuple: Tuple[str, ...],
    before_base: str,
    after_base: str,
    ref_base: str,
    context_model: str,
    fasta_path: str,
    max_matches: int = DEFAULT_MAX_MATCHES
) -> Tuple[Tuple[str, int], ...]:
    """
    Cached version of get_matching_positions.
    Returns tuple for hashability.
    """
    with pysam.FastaFile(fasta_path) as fasta:
        matching = []
        for region in regions_tuple:
            try:
                parts = region.split()
                if len(parts) < 3:
                    logger.warning(f"Invalid region format: {region}")
                    continue
                chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            except (ValueError, IndexError) as e:
                logger.warning(f"Invalid region format: {region}, error: {e}")
                continue
            
            # Ensure we have valid range and avoid boundary issues
            chrom_length = fasta.get_reference_length(chrom)
            start = max(2, start)  # Need at least position 2 for trinucleotide context
            end = min(end, chrom_length - 1)  # Leave room for after base
            
            if start >= end:
                logger.debug(f"Invalid range for {chrom}: start={start}, end={end}")
                continue
            
            for pos in range(start, end):
                try:
                    if is_random_pos_wanted(fasta, chrom, pos, before_base, after_base, ref_base, context_model):
                        matching.append((chrom, pos))
                        if len(matching) >= max_matches:
                            break
                except Exception as e:
                    logger.debug(f"Error checking position {chrom}:{pos}: {e}")
                    continue
                    
            if len(matching) >= max_matches:
                break
                
        return tuple(matching)


def get_matching_positions(
    regions: List[str],
    fasta: pysam.FastaFile,
    before_base: str,
    after_base: str,
    ref_base: str,
    context_model: str,
    max_matches: int = DEFAULT_MAX_MATCHES
) -> List[Tuple[str, int]]:
    """
    Returns a list of all positions in regions that match the given context model.
    Uses cached function for performance.
    """
    # Convert to tuple for caching
    regions_tuple = tuple(regions)
    fasta_path = fasta.filename.decode() if hasattr(fasta.filename, 'decode') else str(fasta.filename)
    
    cached_result = get_matching_positions_cached(
        regions_tuple, before_base, after_base, ref_base, context_model, fasta_path, max_matches
    )
    
    return list(cached_result)


def get_random_mut(
    before_base: str,
    after_base: str,
    ref_base: str,
    regions: List[str],
    fasta: pysam.FastaFile,
    context_model: str,
    model: MutationModel
) -> Tuple[str, int, str, str]:
    """
    Returns a random mutation in a random region that matches the specified criteria.

    Args:
        before_base: base before the position
        after_base: base after the position
        ref_base: reference base
        regions: list of regions in the bed file
        fasta: Fastafile object
        context_model: context model
        model: mutation model

    Returns:
        Tuple of (chromosome, position, reference_base, alternative_base)
        
    Raises:
        ValueError: If no matching context found
    """
    matching_positions = get_matching_positions(
        regions, fasta, before_base, after_base, ref_base, context_model
    )

    if not matching_positions:
        raise ValueError(f"No matching context found for {before_base}-{ref_base}-{after_base} under model {context_model}")

    # Randomly select a position from the matching positions
    random_chr, random_pos = matching_positions[np.random.randint(len(matching_positions))]
    
    # Mutate the base - handle potential None return from mutate
    try:
        alt = model.mutate(ref_base)
        if alt is None:
            raise ValueError(f"Mutation model returned None for base {ref_base}")
    except Exception as e:
        logger.error(f"Error mutating base {ref_base}: {e}")
        raise
    
    return random_chr, random_pos, ref_base, alt


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


def validate_input_files(fasta_path: str, vcf_path: str, bed_path: Optional[str] = None) -> None:
    """Validate that input files exist and are accessible."""
    if not os.path.exists(fasta_path):
        raise FileNotFoundError(f"FASTA file not found: {fasta_path}")
    if not os.path.exists(vcf_path):
        raise FileNotFoundError(f"VCF file not found: {vcf_path}")
    if bed_path and not os.path.exists(bed_path):
        raise FileNotFoundError(f"BED file not found: {bed_path}")


def load_regions(bed_path: Optional[str], fasta: pysam.FastaFile) -> List[str]:
    """Load genomic regions from BED file or generate from FASTA."""
    regions = []
    
    if bed_path is not None:
        try:
            with my_open(bed_path, "r") as f:
                for line in f:
                    line = line.strip()
                    if line and not line.startswith("#"):
                        regions.append(line)
        except Exception as e:
            logger.error(f"Error reading BED file {bed_path}: {e}")
            raise
    else:
        # Generate regions from FASTA references
        for chrom in fasta.references:
            length = fasta.get_reference_length(chrom)
            regions.append(f"{chrom}\t0\t{length}")
    
    logger.info(f"Loaded {len(regions)} regions")
    return regions


def process_vcf_record(
    record, 
    fasta: pysam.FastaFile, 
    regions: List[str], 
    context_model: str, 
    mutation_model: MutationModel
) -> Optional[Dict[str, Any]]:
    """
    Process a single VCF record and return dN/dS statistics.
    
    Returns:
        Dictionary with dN/dS statistics or None if processing failed
    """
    try:
        chrom = record.chrom
        pos = record.pos
        ref_base = record.ref
        alt_base = record.alts[0] if record.alts else None
        
        if alt_base is None:
            logger.warning(f"No alternative allele for {chrom}:{pos}")
            return None
        
        # Get the trinucleotide context
        context_result = get_trinucleotide_context(chrom, pos, fasta)
        if context_result is None:
            logger.warning(f"Invalid trinucleotide context at {chrom}:{pos}")
            return None
            
        before_base, ref_base_context, after_base = context_result
        
        # Find matching random position within the regions
        matching_positions = get_matching_positions(
            regions, fasta, before_base, after_base, ref_base_context, context_model
        )
        
        if not matching_positions:
            logger.warning(f"No matching positions found for {chrom}:{pos}")
            return None
        
        # Select a random position from the matching positions
        random_chr, random_pos = matching_positions[np.random.randint(len(matching_positions))]
        
        # Get codon context
        codon_context = get_trinucleotide_context(random_chr, random_pos, fasta)
        if codon_context is None:
            logger.warning(f"Invalid codon context at {random_chr}:{random_pos}")
            return None
            
        codon = ''.join(codon_context)
        
        # Get the mutated codon - ensure we have valid bases
        try:
            random_alt_base = mutation_model.mutate(ref_base_context)
            if random_alt_base is None:
                logger.warning(f"Mutation model returned None for base {ref_base_context}")
                return None
        except Exception as e:
            logger.warning(f"Error mutating base {ref_base_context}: {e}")
            return None
            
        mutated_codon = codon_context[0] + random_alt_base + codon_context[2]
        
        # Calculate dN/dS statistics for the mutated codon
        return context_dnds(codon, mutated_codon)
        
    except Exception as e:
        logger.error(f"Error processing record {record.chrom}:{record.pos}: {e}")
        return None


def calculate_dnds_ratio(dnds_data: Dict[str, int]) -> float:
    """Calculate final dN/dS ratio from accumulated data."""
    dN = dnds_data["non_synonymous"] / dnds_data["N_sites"] if dnds_data["N_sites"] > 0 else 0
    dS = dnds_data["synonymous"] / dnds_data["S_sites"] if dnds_data["S_sites"] > 0 else 0
    return dN / dS if dS > 0 else float('inf')


def context(
    fasta: str, 
    vcf: str, 
    bed: str = None, 
    model: str = "random", 
    alpha: float = None, 
    beta: float = None, 
    gamma: float = None, 
    pi_a: float = None,
    pi_c: float = None,
    pi_g: float = None,
    pi_t: float = None,
    context_model: str = "codon"
) -> float:
    """
    Run the simulation and calculate the dN/dS ratio.
    
    Args:
        fasta: Path to FASTA file
        vcf: Path to VCF file
        bed: Path to BED file (optional)
        model: Mutation model type
        alpha, beta, gamma: Mutation model parameters
        pi_a, pi_c, pi_g, pi_t: Base frequency parameters
        context_model: Context model type
        
    Returns:
        dN/dS ratio
        
    Raises:
        FileNotFoundError: If input files don't exist
        ValueError: If invalid parameters provided
    """
    # Validate input files
    validate_input_files(fasta, vcf, bed)
    
    # Open files with proper resource management
    with pysam.FastaFile(fasta) as fasta_file, \
         pysam.VariantFile(vcf) as vcf_file:
        
        # Load regions
        regions = load_regions(bed, fasta_file)
        
        # Create mutation model
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
        
        # Initialize data collection
        dnds_data = defaultdict(int)
        processed_records = 0
        
        # Process VCF records
        for record in vcf_file:
            logger.info(f"Processing record: {record.chrom}:{record.pos} {record.ref} -> {record.alts}")
            
            dnds_stats = process_vcf_record(record, fasta_file, regions, context_model, mutation_model)
            
            if dnds_stats is not None:
                dnds_data["N_sites"] += dnds_stats["N_sites"]
                dnds_data["S_sites"] += dnds_stats["S_sites"]
                if dnds_stats["mutation_type"] == "synonymous":
                    dnds_data["synonymous"] += 1
                else:
                    dnds_data["non_synonymous"] += 1
                processed_records += 1
        
        logger.info(f"Successfully processed {processed_records} records")
        
        # Calculate and return dN/dS ratio
        return calculate_dnds_ratio(dnds_data)

if __name__ == "__main__":
    pass