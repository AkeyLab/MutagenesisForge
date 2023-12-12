# Assuming that bed file and fasta files are provided

# Rule All
rule all:
    input:
        "stats.txt"

# Run simulations
rule sim:
    params:
        vcf=config['vcf'],
        bed=config['bed'],
        fasta=config['fasta'],
        sims=config['sims']
    output:
        expand("vep_output_{index}.txt", index=range(1, lambda wildcards: params.num_output_files + 1))
        # "vep.output"
    shell:
        "src/sim.py --fasta {params.fasta} --input_bed {params.bed} --input_mut {params.vcf}, --sim_count {params.sims}"

# Run stats
rule stats:
    input:
        lambda wildcards: expand("vep_output_{index}.txt", index=range(int(wildcards.arg) + 1))
        # "vep.output"
    params:
        stat=config['stat']
    output:
        "stats.txt"
    shell:
        "src/stats.py --stat {params.stat}"

