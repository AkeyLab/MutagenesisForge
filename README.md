# MutagenesisForge
CLI and python package of tools written to work with randomized mutation simulations to generate variant effects and calculate simulated transition-transversion ratios

## installation
```bash
pip install git+https://https://github.com/AkeyLab/MutagenesisForge
```

## Command Line Interface

### vcf-construction
Given input files, a VCF (variant call format) file is generated of randomized mutations with matching trinucleotide context, within the bed regions of the input fasta file.

#### Inputs
- vcf:
vcf file path of original mutations
- bed:
bed file path to be used 
- fasta:
fasta file path used 
- sims:
number of vcf files generated (default = 1)
- tstv:
transition-transversion ratio used for mutation generation (default = 2)
- output:
prefix for output file (default = 'output.vcf')
- vep-call:
boolean which determines if Variant Effect Predictor (VEP) software is to be run on each of the created VCF files.

Example usage:
```bash
MutagenesisForge vcf-construction --vcf ex.vcf --bed ex.bed --fasta ex.fa --tstv 2.5 --sims 40
```
returns vcf file of randomizd mutations generated from random mutation of trinucleotide context of input vcf as found within bed file regions in fasta file

### exhaust
Given an input fasta file, a transition-transversion rate is calculated based off of all possible permutations of possible codons.

#### Inputs
- fasta:
fasta file path used
- by-read:
tstv calculation is done by a mean of each read's tstv (default = False)

Example usage:
```bash
MutagenesisForge exhaust --fasta ex.fa
```
returns tstv value of each permutation of all codon positions of fasta file ex.fa

### tstv-test
A transition-transversion ratio is calculated from an input VCF file

#### Inputs
- vcf:
vcf file path used

Example usage:
```bash
MutagenesisForge tstv-test ex.vcf
```
returns the transition-transversion ratio of file ex.vcf

## Python Package Documentation
TODO: To be written...
