# from .utils import *
import click
#from scipy import stats
from .random_vcf import vcf_constr
from .stat_utils import synoymous_detected, dNds, total_dNds
from .transition_transversion import tstv
from .exhaustive import exhaustive

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
    '--model',
    type=click.Choice(['random', 'K2P', 'K3P']),
    default = 'random',
    help='evolutionary model'
)
@click.option(
    '--output',
    default = 'output.vcf',
    help='output name'
)
@click.option(
    '--sims',
    default = 1,
    help='number of simulations'
)
@click.option(
    '--vep_call',
    default = False,
    help='run vep call'
)
@click.option(
    '--alpha',
    default = None,
    help='alpha parameter' 
)
@click.option(
    '--beta',
    default = None,
    help='beta parameter'
)
@click.option(
    '--gamma',
    default = None,
    help='gamma parameter'
)
# context model
def context(vcf, bed, fasta, output, sims, model, alpha, beta, gamma, vep_call):
    """
    Given a bed file, mutation file, fasta file, output
    file, transition-transversion ratio, and number of simulations,
    create a vcf file of random mutations

    Parameters:
        bed (str): path to bed file
        vcf (str): path to vcf file
        fasta (str): path to fasta file
        output (str): output file name
        tstv (float): transition-transversion ratio
        sims (int): number of simulations
        vep_call (bool): run vep call

    Returns:
        None (creates a vcf file of random mutations)
    """
    click.echo('vcf construction started')
    vcf_constr(bed, vcf, fasta, output, sims, vep_call, model, alpha, beta, gamma)
    print("vcf construction complete")

    # returning dN/dS values if vep_call is True
    if vep_call:
        veps = []
        for i in range(sims):
            vep_path = output + str(i) + '.vep'
            veps.append(vep_path)
        #dnds differences can be calculated here
        print(f'dN/dS for empirical vcf: {dNds(vcf)}')
        print(f'dN/dS for total simulation data: {total_dNds(veps)}')
        print(f'dN/dS for each simulation: {[dNds(vep) for vep in veps]}')


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
def exhaust(fasta, by_read = False):
    """
    Given a fasta file, calculate the dN/dS ratio using exhaustive method 
    where each permutation of the codon is tested

    Parameters:
        fasta (str): path to fasta file
        by_read (bool): calculate dN/dS by gene
    
    Returns:
        None (prints dN/dS ratio)
    """
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
def tstv(vcf):
    """
    TODO: may be removed in futurep
    """
    click.echo('Transition-transversion model started')
    click.echo(f"tstv = {tstv(vcf)}")
    print("Transition-transversion model finished")

if __name__ == '__main__':
    cli()
