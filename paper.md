---
title: 'MutagenesisForge: A framework for modeling codon-level mutational biases and dN/dS selection'
tags:
  - bioinformatics
  - evolutionary genomics
  - codon models
  - dN/dS
  - mutation simulation
  - VEP
authors:
  - name: Cooper Koers
    orcid: 0009-0008-5214-917X
    corresponding: true
    affiliation: 1 # (Multiple affiliations must be quoted)
  - name: Rob F. Bierman
    orcid: 0000-0001-8513-7425
    affiliation: 1
  - name: Huixin Xu
    orcid: 0000-0003-4098-9901
    affiliation: 1
  - name: Joshua M. Akey
    orcid: 0000-0002-4411-1330
    affiliation: 1

affiliations:
 - name: The Lewis-Sigler Institute for Integrative Genomics, Princeton University, Princeton, NJ 08540, USA.
   index: 1
date: 17 Aprul 2025
bibliography: paper.bib
---
# Introduction and Statement of Need
In molecular evolution, point mutations within coding regions are typically classified as either synonymous – those that do not change amino acid sequence – or synonymous, which do. The ratio of nonsynonymous to synonymous substitutions (dN/dS) has become a cornerstone analysis for evaluating selective pressure trends in comparative genomics. A ratio of 1 suggests neutral evolution, <1 implies purifying selection, and >1 imparts positive selection. 

Testing of broad evolutionary pressure utilizing dN/dS ratio analysis
is a practiced assay within bioinformatics, specifically within comparative evolutionary genomics. The utility of such ratios lies in its ability for comparison. Relationships between dN/dS ratios of different datasets offer a means of comparing relative selection strengths. This offers comparative genomic analyses to be made between dataset dN/dS ratios. Despite dN/dS ratios providing a powerful metric for detecting selection, interpretation in isolation can lead to deeply misleading results. While multiple applications exist to perform the analysis, they do not offer comparative models for individual inputs. Without context, the value of this ratio may lack definitive biological meaning.

Despite the availability of tools to calculate dN/dS ratios, popular methods of dN/dS calculation lack a modular framework to produce simulated datasets tailored to empirical data. In response, we present MutagenesisForge: a command-line tool designed to produce simulated results to aid in comparative genomic analysis. MutagenesisForge allows for the construction of input data-built codon-specific null models of evolution designed to offer data complementary to dN/dS values produced from popular calculator tools. With MutagenesisForge, simulation functionality is brought to the user through one command line interface.

# Design
MutagenesisForge offers two distinct methods to present simulated data to users: Exhaustive and Contextual, both of which allow comparison between empirical dN/dS ratios to null model data developed from user-inputted data and evolutionary model.

## Exhaustive Model
The exhaustive method systematically simulates all possible single-nucleotide substitutions across an input reference FASTA file. Beginning at the first found start codon, the algorithm consists of a sliding window analysis incrementing by codon, where each position has each substitution performed and categorized as either synonymous or nonsynonymous. The resulting counts of mutation type are then divided to return the simulated dN/dS ratio. This method is compatible with files with FASTA format and offers user selection for several evolutionary substitution models. Users are provided with the optional flag to restrict the exhaustive search to genomic regions of interest.

## Context Model
The context simulation model draws from the context of Variant Call Format (VCF) nucleotides by forming codons from bases flanking the nucleotides of interest from a reference FASTA file. This codon is then stochastically mapped back to regions of interest within reference genome files. After mapping, the nucleotide of interest is then mutated to a base at a certain probability given the evolutionary model selected by the user. Each nucleotide from the original VCF will have its simulated coordinate and mutation stored in an output VCF file for downstream analysis of dN/dS calculation.

### Variant Effect Predictor Integration
Variant Effect Predictor Integration 
A feature of the context model is the integration of Variant Effect Predictor (VEP) from ENSMBL. Simulated VCF files of  generated from the Context model may automatically be passed to VEP for annotation, providing users with helpful to-date annotations of supported species.

# Applications in Research
MutagenesisForge was used in Xu et al. 2025, where it served in analysis of new somatic mutation files to compare between statistical methods of searching for de novo mutations found in other publications. MutagenesisForge-produced data served as a null model for comparison of statistically-derived mutation datasets. By comparing the mutation spectrum of actual mutation-finding models, the interface allowed for reinforced statistical assessment.