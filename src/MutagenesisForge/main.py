# from .utils import *
import click
#from scipy import stats
from .context import context as dnds_context
from .stat_utils import synoymous_detected, dNds, total_dNds
from .transition_transversion import tstv
from .exhaustive import exhaustive

"""
this file is the main entry point for the MutagenesisForge package
"""

# group cli test options
@click.group()
def cli():
    pass

# click command for context method
@cli.command()
@click.option(
    '--vcf',
    prompt='path to reference vcf file',
    help='path to reference vcf file',
    required=True
)
@click.option(
    '--bed',
    '-b',
    prompt='path to bed file',
    help='path to bed file',
    default=None
)
@click.option(
    '--fasta',
    '-f',
    prompt='path to fasta file',
    help='path to fasta file',
    required=True
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
@click.option(
    '--pi_a',
    default = None,
    help='pi_a parameter for HKY85 and F81 models'
)
@click.option(
    '--pi_c',
    default = None,
    help='pi_c parameter for HKY85 and F81 models'
)
@click.option(
    '--pi_g',
    default = None,
    help='pi_g parameter for HKY85 and F81 models'
)
@click.option(
    '--pi_t',
    default = None,
    help='pi_t parameter for HKY85 and F81 models'
)
@click.option(
    '--context_model',
    type = click.Choice(['blind', 'ra', 'ra_ba', 'ra_aa', 'codon']),
    default = 'codon',
    help='context model used for mutation'
)
@click.option(
    '--model',
    '-m',
    type=click.Choice(['random', 'JC69', 'K2P', 'F81', 'HKY85']),
    default='random',
    help='evolutionary model to use for context calculation'
)
@click.option(
    '--sims',
    '-s',
    default=1,
    help='Number of simulations to run'
)
# verbose flag
@click.option(
    '--verbose',
    '-v',
    is_flag=True,
    help='Print verbose output'
)
# context model
def context(
    vcf,
    bed,
    fasta,
    model,
    output,
    alpha,
    beta,
    gamma,
    pi_a,
    pi_c,
    pi_g,
    pi_t,
    context_model,
    sims,
    verbose):
    """
    Return a dN/dS ratio given the context data provided

    Parameters:
        vcf (str): path to vcf file
        bed (str): path to bed file
        fasta (str): path to fasta file
        alpha (float): alpha parameter
        beta (float): beta parameter
        gamma (float): gamma parameter
        pi_a (float): pi_a parameter
        pi_c (float): pi_c parameter
        pi_g (float): pi_g parameter
        pi_t (float): pi_t parameter
        context_model (str): context model
        model (str): evolutionary model to use for context calculation
        sims (int): number of simulations to run
        verbose (bool): if True, print verbose output

    Returns:
        None (prints dN/dS ratio)
    """
    if verbose:
        click.echo('Verbose mode enabled')
        click.echo('Context model started')
        click.echo(f"Calculating dN/dS ratio using {context_model} model with {model} evolutionary model")
    for i in range(sims):
        if verbose:
            click.echo(f"Simulation {i+1}/{sims}")
        # Call the dnds_context function with the provided parameters
        dnds = dnds_context(
            vcf=vcf, 
            bed=bed, 
            fasta=fasta, 
            alpha=alpha, 
            beta=beta, 
            gamma=gamma, 
            pi_a=pi_a, 
            pi_c=pi_c, 
            pi_g=pi_g, 
            pi_t=pi_t, 
            context_model=context_model,
            model=model
        )
        if verbose:
            click.echo(f"dN/dS ratio for simulation {i+1}: {dnds}")
        else:
            click.echo(f"{dnds}")


# click command for exhaustive method

"""
(fasta: str, bed: str | None = None, by_read: bool = False, model: str = "random", alpha: str = None, beta: str = None, gamma: str = None, pi_a: str = None, pi_c: str = None, pi_g: str = None, pi_t: str = None) -> (floating[Any] | float)
"""
@cli.command()
@click.option(
    '--fasta',
    '-f',
    prompt='Path to fasta file',
    help='Path to fasta file',
)
@click.option(
    '--bed',
    '-b',
    prompt='Path to bed file',
    default=None,
    help='Path to bed file',
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
@click.option(
    '--pi_a',
    default = None,
    help='pi_a parameter for HKY85 and F81 models'
)
@click.option(
    '--pi_c',
    default = None,
    help='pi_c parameter for HKY85 and F81 models'
)
@click.option(
    '--pi_g',
    default = None,
    help='pi_g parameter for HKY85 and F81 models'
)
@click.option(
    '--pi_t',
    default = None,
    help='pi_t parameter for HKY85 and F81 models'
)
@click.option(
    '--model',
    '-m',
    type=click.Choice(['random', 'JC69', 'K2P', 'F81', 'HKY85']),
    default='random',
    help='evolutionary model to use for context calculation'
)
@click.option(
    '--sims',
    '-s',
    default=1,
    help='Number of simulations to run'
)
# verbose flag
@click.option(
    '--verbose',
    '-v',
    is_flag=True,
    help='Print verbose output'
)
# flag to calculate dN/dS by gene
@click.option(
'--by-read',
is_flag=True,
help='Calculate dN/dS by gene'
)
def exhaust(fasta, bed, alpha, beta, gamma, pi_a, pi_c, pi_g, pi_t, model, sims, verbose, by_read):
    """
    Given a fasta file, calculate the dN/dS ratio using exhaustive method 
    where each permutation of the codon is tested

    Parameters:
        fasta (str): path to fasta file
        bed (str): path to bed file
        alpha (float): alpha parameter
        beta (float): beta parameter
        gamma (float): gamma parameter
        pi_a (float): pi_a parameter for HKY85 and F81 models
        pi_c (float): pi_c parameter for HKY85 and F81 models
        pi_g (float): pi_g parameter for HKY85 and F81 models
        pi_t (float): pi_t parameter for HKY85 and F81 models
        model (str): evolutionary model to use for context calculation
        sims (int): number of simulations to run
        verbose (bool): if True, print verbose output
        by_read (bool): if True, calculate dN/dS by gene
    Returns:
        None (prints dN/dS ratio)
    """
    if by_read:
        if verbose:
            click.echo('Verbose mode enabled')
            click.echo('Exhaustive model started')
            click.echo(f"Calculating dN/dS ratio using {model} evolutionary model")
        dnds = exhaustive(fasta, bed, alpha, beta, gamma, pi_a, pi_c, pi_g, pi_t, model, sims, verbose, by_read=True)
        if verbose:
            click.echo(f"dN/dS ratio for each gene: {dnds}")
        else:
            click.echo(f"{dnds}")
    else:
        if verbose:
            click.echo('Verbose mode enabled')
            click.echo('Exhaustive model started')
            click.echo(f"Calculating dN/dS ratio using {model} evolutionary model")
        dnds = exhaustive(fasta, bed, alpha, beta, gamma, pi_a, pi_c, pi_g, pi_t, model, sims, verbose, by_read=False)
        if verbose:
            click.echo(f"dN/dS = {dnds}")
        else:
            click.echo(f"{dnds}")

if __name__ == '__main__':
    cli()
