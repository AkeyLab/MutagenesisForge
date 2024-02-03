from .sim import sim
from .stat_utils import synoymous_detected, dNds, total_dNds
from .exhaustive import exhaustive
from .transition_transversion import tstv
from .utils import load_parameter_from_yaml

import os
import click
import numpy as np
from scipy import stats


        
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

# group cli test options
@click.group()
def cli():
    pass

# click command for sim method
@cli.command()
def sim_method():
    print("Sim model")

# click command for exhaustive method
@cli.command()
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
def exhaustive_method(fasta, by_read = False):
    if by_read:
        click.echo("Exhaustive model ratio of each gene")
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
