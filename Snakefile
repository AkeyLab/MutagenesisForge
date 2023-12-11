# Assuming that bed file and fasta files are provided

# Rule All
rule all:
    input:
        "stats.txt"

# Run simulations
rule sim:
    input:
        "ex.vcf",
        "ex.bed",
        "ex.fa"
    output:
        expand("vep_output_{index}.txt", index=range(1, lambda wildcards: params.num_output_files + 1))
        # "vep.output"
    shell:
        "src/sim.py --fasta ex.fa --input_bed ex.bed --input_mut ex.vcf, --sim_num 2"

# Run stats
rule stats:
    input:
        lambda wildcards: expand("vep_output_{index}.txt", index=range(int(wildcards.arg) + 1))
        # "vep.output"
    output:
        "stats.txt"
    shell:
        "touch stats.txt"

