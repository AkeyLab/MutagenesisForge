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
        output_files=dynamic(lambda wildcards: ["vep_output_{}.txt".format(i) for i in range(1, int(wildcards.params.sims) + 1)])
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
        "stats_{params.stat}.txt"
    shell:
        "src/stats.py --stat {params.stat}"

