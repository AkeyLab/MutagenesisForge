from .sim import sim
from .stat_utils import synoymous_detected, dNds, total_dNds
from .exhaustive import exhaustive
from .transition_transversion import tstv
from .utils import load_parameter_from_yaml

import os
import click
import numpy as np
from scipy import stats

from vep_sim.mut_vcf import vcf_constr


# group cli test options
@click.group()
def cli():
    pass

# click command for vcf construction debugging
@cli.command()
@click.option(
    '--vcf',
    prompt='path to reference vcf file',
    help='path to reference vcf file'
)
@click.option(
    '--bed',
    prompt='path to bed file',
    help='path to bed file'
)
@click.option(
    '--fasta',
    prompt='path to fasta file',
    help='path to fasta file'
)
@click.option(
    '--tstv',
    default = 2.0,
    help='transition-transversion ratio'
)
@click.option(
    '--output',
    default = 'output.vcf',
    help='output name'
)
def vcf_construction(vcf, bed, fasta, output, tstv = 2.0):
    click.echo('vcf construction started')
    if tstv != 2.0:
        vcf_constr(vcf, bed, fasta, output, tstv)
    if tstv == 2.0:
        vcf_constr(vcf, bed, fasta, output) 
    print("vcf construction complete")

# click command for vep debugging
@cli.command()
@click.option(
    '--vcf',
    prompt='path to reference vcf file',
    help='path to reference vcf file'
)
# click command for sim method
@cli.command()
@click.option(
    '--vcf',
    prompt='path to reference vcf file',
    help='path to reference vcf file'
)
@click.option(
    '--bed',
    prompt='path to bed file',
    help='path to bed file'
)
@click.option(
    '--fasta',
    prompt='path to fasta file',
    help='path to fasta file'
)
@click.option(
    '--sim-num',
    prompt='number of simulations',
    default=10,
    help='number of simulations'
)
@click.option(
    '--output',
    default = 'output.vep',
    help='output name'
)
@click.option(
    '--tstv',
    default = 2.0,
    help='transition-transversion ratio'
)
def sim_method(vcf, bed, fasta, sim_num, output, tstv = 2.0):
    click.echo('Sim model started')
    if tstv != 2.0:
        sim(vcf, bed, fasta, sim_num, output, tstv)
    if tstv == 2.0:
        sim(vcf, bed, fasta, sim_num, output) 
    print("Sim model complete")

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
