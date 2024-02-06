from .sim import sim
from .stat_utils import synoymous_detected, dNds, total_dNds
from .exhaustive import exhaustive
from .transition_transversion import tstv
from .utils import load_parameter_from_yaml

import os
import click
import numpy as np
from scipy import stats
import tempfile


        
# apply this function to the directory of vep output files for statistical analysis
# only used in sim method
#def find_files_with_prefix(prefix, directory = load_parameter_from_yaml('parameters.yaml', 'output_dir')):
#    matching_files = []
#    for filename in os.listdir(directory):
#        if filename.startswith(prefix):
#            matching_files.append(filename)
#    return matching_files
        

# access parameters from yaml file for sim method
#emp_vep_path = load_parameter_from_yaml('parameters.yaml', 'emp_vep_path')
#sim_vep_path = load_parameter_from_yaml('parameters.yaml', 'sim_vep_path')
#output_dir = load_parameter_from_yaml('parameters.yaml', 'output_dir')

# potential functionally for slurm job submission
#def submit_job(file):
#    os.system(f'sbatch {file}')
#    return

# example of sbatch file submission without .sh file
#sbatch <<EOF
## Slurm directives and commands
##SBATCH --job-name=my_job
##SBATCH --output=output.txt
##SBATCH --error=error.txt
##SBATCH --partition=compute
##SBATCH --nodes=1
##SBATCH --ntasks=4
##SBATCH --cpus-per-task=1
##SBATCH --mem=4G
##SBATCH --time=1:00:00

# Commands to run
#your_command arg1 arg2
#EOF

def slurm_submit(command):
    with tempfile.TemporaryDirectory() as temp_dir:
        f = os.path.join(temp_dir, "slurm_test.sh")
        with open("slurm_test.sh", "w") as f:
            f.write("#!/bin/bash\n")
            f.write("#SBATCH --job-name=test\n")
            f.write("#SBATCH --output=output.txt\n")
            f.write("#SBATCH --error=error.txt\n")
            f.write("#SBATCH --partition=compute\n")
            f.write("#SBATCH --nodes=1\n")
            f.write("#SBATCH --ntasks=4\n")
            f.write("#SBATCH --cpus-per-task=1\n")
            f.write("#SBATCH --mem=4G\n")
            f.write("#SBATCH --time=1:00:00\n")
            f.write("\n")
            f.write("your_command arg1 arg2\n")
            os.system(f'sbatch {f}')


# group cli test options
@click.group()
@click.option(
    '--slurm',
    is_flag=True,
    help='Submit job to slurm'
)
def cli():
    pass

# click command for sim method
@cli.command()
def sim_method():
    print("Sim model")

# click command for exhaustive method
@cli.command()

# click option for fasta file
@click.option(
    '--fasta',
    prompt='Path to fasta file',
    help='Path to fasta file',
)

# flag to calculate dN/dS by gene
@click.option(
'--by-read',
is_flag=True,
help='Calculate dN/dS by gene'
)

# method for exhaustive method
def exhaustive_method(fasta, by_read = False, slurm = False):
    if by_read:
        # if read flag as additional conditional perhaps
        click.echo("Exhaustive model mean ratio of each read")
        dnds = exhaustive(fasta, by_read=True)
        click.echo(f"dN/dS = {dnds}")
        click.echo("Exhaustive model complete")
    else:
        click.echo("Exhaustive model ratio of entire file")
        dnds = exhaustive(fasta)
        click.echo(f"dN/dS = {dnds}")
        click.echo("Exhaustive model complete")

# click command for tstv method
@cli.command()
@click.option(
    '--vcf',
    prompt='Path to vcf file',
    help='Path to vcf file',
)
def tstv_test(vcf):
    click.echo('Transition-transversion model started')
    click.echo(f"tstv = {tstv(vcf)}")
    print("Transition-transversion model finished")

if __name__ == '__main__':
    cli()
