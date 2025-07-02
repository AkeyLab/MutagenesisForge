

<p align="center">
  <img src="logo.png" width="300" alt="MutagenesisForge Logo" />
</p>

------------


# MutagenesisForge - A framework for modeling codon-level mutational biases and d~N~/d~S~ selection

## Overview
**MutagenesisForge** is a codon-level simulation software designed to calculate d~N~/d~S~ null models built on user-inputed data for comparative genomic analyses. It allows for the construction of input data-built codon-specific null models of evolution designed to offer data complementary to d~N~/d~S~ values produced from popular calculator tools. With MutagenesisForge, simulation functionality is brought to the user through one command line interface.

Simulation Methods:

- Exhaustive: 
- Contextual: 


## Installation
```bash
pip install git+https://https://github.com/AkeyLab/MutagenesisForge
```
## Exhaustive Codon Model Example
```
MutagenesisForge exhaust -fasta
```


## Context Codon Model
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
