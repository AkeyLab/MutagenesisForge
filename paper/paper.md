---
title: 'MutagenesisForge: A framework for modeling codon-level mutational biases and dN/dS selection'
tags:
  - bioinformatics
  - Python
  - evolutionary genomics
  - codon models
  - dN/dS
  - mutation simulation
  - comparative genomics
authors:
  - name: Cooper Koers
    orcid: 0009-0008-5214-917X
    corresponding: true
    affiliation: 1 
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
date: 13 July 2025
bibliography: paper.bib
---

# Introduction and Statement of Need

Point mutations that change the DNA sequence in protein-coding regions of the genome are typically classified as synonymous, nonsynonymous (or missense), or protein-truncating (or nonsense) depending on their functional consequences [@nei1986]. Synonymous mutations do not change the amino acid sequence of a protein whereas nonsynonymous do result in amino acid changes. Protein-truncating mutations create a stop codon that results in translation of the protein to prematurely end. We focus on synonymous and nonsynonymous mutations, which are far more abundant in catalogs of natural protein-coding sequence variants and because they have become a cornerstone in evolutionary analyses. Specifically, the ratio of nonsynonymous to synonymous substitutions (d~N~/d~S~) between two or more species is a fundamentally important measure that evaluates the selective pressures that a protein has experienced [@williams2020]. A d~N~/d~S~ = 1 suggests neutral evolution, d~N~/d~S~ < 1 implies purifying selection, and d~N~/d~S~ > 1 is consistent with the effects of positive selection [@kryazhimskiy2008]. Comparing d~N~/d~S~ ratios across protein-coding genes and species provides considerable insights into the types of proteins and biological processes that have been adaptive in the evolutionary history of distinct species [@bustamante2005; @nielsen2005]. Comparisons of d~N~/d~S~ ratios across datasets enable the evaluation of relative selection strengths. This offers comparative genomic analyses to be made between dataset d~N~/d~S~ ratios. Despite d~N~/d~S~ ratios providing a powerful metric for detecting selection, interpretation in isolation can lead to deeply misleading results [@hughes2007]. Additionally important is the substitution context of datasets. Differing models of substitution will yield often incomparable ratios [@goldman1994]. While multiple applications exist to perform the analysis, they do not offer comparative models for individual inputs or the ability to account for evolutionary pretenses. Without context, the value of this ratio may lack definitive biological meaning.

Despite the availability of tools to calculate d~N~/d~S~ ratios, popular methods of d~N~/d~S~ calculation lack a modular framework to produce simulated datasets tailored to empirical data [@pond2005; @yang2007; @zhang2006]. In response, we present `MutagenesisForge`: a command-line tool designed to produce simulated results to aid in comparative genomic analysis. `MutagenesisForge` allows for the construction of input data-built codon-specific null models of evolution designed to offer data complementary to d~N~/d~S~ values produced from popular calculator tools. With `MutagenesisForge`, simulation functionality is brought to the user directly through either a command line interface (CLI) or as an installable python package for programmatic use.

# Design

## Evolutionary Model-Dependant Substitution for Codon-level Mutational Biases

`MutagenesisForge` offers a platform to simulate nucleotide-level substitutions according to user-inputted evolutionary models. The models are shared with the exhaustive and context methods, allowing for continuity across the framework. Given an input nucleotide and substitution model with given parameters, a corresponding nucleotide is returned following the probabilities of the substitution matrix.

![**Matrix representations of supported evolutionary substitution models.**: MutagenesisForge supports the Kimura 2-parameter (K2P), Kimura 3-Parameter (K3P), Jukes-Cantor (JC69), Felenstein 4-parameter (F81), and Hasegawa, Kishino, and Yano (HKY85) substitution models of evolution to simulate single-nucleotide substitution. Rows and columns of each matrix represent the nucleotides thymine, cytosine, adenine, and guanine.](fig1.png) Matrix representations of supported evolutionary substitution models [@felsenstein1981; @hasegawa1985; @jukes1969; @kimura1980; @kimura1981].


## d~N~/d~S~ Simulation Methods

MutagenesisForge offers two distinct methods to present simulated data to users: `Exhaustive` and `Contextual`, both of which allow comparison between empirical d~N~/d~S~ ratios to null model data developed from user-inputted data and evolutionary model.

### Exhaustive Model

The exhaustive method systematically simulates all possible single-nucleotide substitutions across an input reference FASTA file. Starting from the first identified start codon, the algorithm consists of a sliding window analysis incrementing by codon, where each position has each substitution performed and categorized as either synonymous or nonsynonymous. The resulting counts of mutation type and mutation type site are used as inputs to calculate the d~N~/d~S~ ratio. This method accepts FASTA files and supports user-defined evolutionary substitution models. Users are provided with the optional flag to restrict the exhaustive search to genomic regions of interest.

![**Generating d~N~/d~S~ values using `MutagenesisForge` Exhaustive Model According to K2P Substitution Model.**: The logic displayed represents the manner `MutagenesisForge` calculates the number of missense and synonymous per codon in accordance with the respective substitution model.](fig2.png)

### Context Model

The context model simulates mutations based on the nucleotide context surrounding variants in a VCF file by forming codons from bases flanking the nucleotides of interest from a reference FASTA file [@aggarwala2016]. This codon is then stochastically mapped back to regions of interest within reference genome files. After mapping, the nucleotide of interest is then mutated to a base at a certain probability given the evolutionary model selected by the user. Analysis between the mapped codon and mutated mapped codon then allows for the calculation of null model d~N~/d~S~ from the lens of the nucleotide context of the input file.

![**Finding Contextually-Valid Codons and Performing Evolutionary-Aware Mutation Using MutagenesisForge Context.**: Given the input matched context of the base and upstream base, a stochastic match is found within the coordinates of interest. From this found codon, a substitution is made according to the input substitution model. Data from each substitution is then used to calculate the d~N~/d~S~ value of the simulation to return a context-aware simulated value.](fig3.png)

## Applications in Research

MutagenesisForge was used in Xu et al. 2025, where it served in analysis of new somatic mutation files to compare between statistical methods of searching for de novo mutations found in other publications. MutagenesisForge-produced data served as a null model for comparison of statistically-derived mutation datasets. By comparing the mutation spectrum of actual mutation-finding models, the interface allowed for reinforced statistical assessment.

# Acknowledgements

We acknowledge the support and intellectual input of members of the Akey Lab. We also thank the Lewis-Sigler Institute for Integrative Genomics for harboring a collaborate and welcoming research enviornment. This work was not supported by external funding.

# References
