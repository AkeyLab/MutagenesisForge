from .sim import sim
from .stat_utils import synoymous_detected, dNds, total_dNds
from .exhaustive import exhaustive
import os
import click
import numpy as np
from scipy import stats
from .utils import load_parameter_from_yaml

        
# apply this function to the directory of vep output files for statistical analysis
def find_files_with_prefix(prefix, directory = load_parameter_from_yaml('parameters.yaml', 'output_dir')):
    matching_files = []
    for filename in os.listdir(directory):
        if filename.startswith(prefix):
            matching_files.append(filename)
    return matching_files
        

# access parameters from yaml file
emp_vep_path = load_parameter_from_yaml('parameters.yaml', 'emp_vep_path')
sim_vep_path = load_parameter_from_yaml('parameters.yaml', 'sim_vep_path')
output_dir = load_parameter_from_yaml('parameters.yaml', 'output_dir')

@click.group()
def cli():
    pass

@cli.command()
def sim_method():
    pass

@cli.command()
@click.argument(
    'fasta_path'
)
def exhaustive_method(fasta_path):
    exhaustive.exhaustive(fasta_path, output_dir)
    print("Exhaustive model")

@cli.command()
def transition_transversion():
    print("Transition transversion model")

if __name__ == '__main__':
    cli()