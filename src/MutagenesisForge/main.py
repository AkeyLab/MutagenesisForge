import click
from .context2 import context as dnds_context
from .exhaustive import exhaustive as dnds_exhaustive

"""
this file is the main entry point for the MutagenesisForge package
"""

"""
def context(    
    fasta: str, 
    vcf: str,
    model: str = "random", 
    alpha: float = 2.0, 
    beta: float = 1.0, 
    gamma: float = 1.0, 
    pi_a: float = 0.3,
    pi_c: float = 0.2,
    pi_g: float = 0.2,
    pi_t: float = 0.3,
    omega: float = 0.5)

    """

# group cli test options
@click.group()
def cli():
    pass

# click command for context method
@cli.command()
@click.option(
    '--fasta',
    '-f',
    help='path to fasta file',
    required=True
)
@click.option(
    '--vcf',
    '-v',
    help='path to reference vcf file',
    required=True
)
@click.option(
    '--alpha',
    default = 2.0,
    help='alpha parameter',
    required=False
)
@click.option(
    '--beta',
    default = 1.0,
    help='beta parameter',
    required=False
)
@click.option(
    '--gamma',
    default = 1.0,
    help='gamma parameter',
    required=False
)
@click.option(
    '--pi_a',
    default = 0.3,
    help='pi_a parameter for HKY85 and F81 models',
    required=False
)
@click.option(
    '--pi_c',
    default = 0.2,
    help='pi_c parameter for HKY85 and F81 models',
    required=False
)
@click.option(
    '--pi_g',
    default = 0.2,
    help='pi_g parameter for HKY85 and F81 models',
    required=False
)
@click.option(
    '--pi_t',
    default = None,
    help='pi_t parameter for HKY85 and F81 models',
    required=False
)
@click.option(
    '--omega',
    default = 0.5,
    help='omega parameter for K3P model',
    required=False
)
@click.option(
    '--model',
    '-m',
    type=click.Choice(['random', 'JC69', 'K2P', 'F81', 'HKY85', 'K3P']),
    default='JC69',
    help='evolutionary model to use for context calculation',
    required=True
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
    '-V',
    is_flag=True,
    help='Print verbose output'
)
# context model
def context(
    vcf,
    fasta,
    model,
    alpha,
    beta,
    gamma,
    pi_a,
    pi_c,
    pi_g,
    pi_t,
    omega,
    sims,
    verbose):
    """
    Return a dN/dS ratio given the context data provided
    """
    ratios = []
    if verbose:
        click.echo('Verbose mode enabled')
        click.echo('Context model started')
        click.echo(f"Calculating dN/dS ratio using evolutionary model")
    for i in range(sims):
        if verbose:
            click.echo(f"Simulation {i+1}/{sims}")
        # Call the dnds_context function with the provided parameters
        dnds = dnds_context(
            vcf=vcf,  
            fasta=fasta, 
            alpha=alpha, 
            beta=beta, 
            gamma=gamma, 
            pi_a=pi_a, 
            pi_c=pi_c, 
            pi_g=pi_g, 
            pi_t=pi_t, 
            omega=omega,
            model=model
        )
        ratios.append(dnds)
        if verbose:
            click.echo(f"dN/dS ratio for simulation {i+1}: {dnds}")
        else:
            click.echo(f"{dnds}")
    if verbose:
        click.echo(f"All dN/dS ratios: {ratios}")
    else:
        click.echo(f"{ratios}")

# click command for exhaustive method

"""
def exhaustive(fasta: str, 
               bed: Optional[str] = None, 
               by_read: bool = False, 
               model:str = "random", 
               alpha:str = 2.0, 
               beta:str = 1.0, 
               gamma:str = 1.0,
               pi_a:str = 0.3,
               pi_c:str = 0.2,
               pi_g:str = 0.2,
               pi_t:str = 0.3,
               omega:str = 0.5) -> float:
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
    default=None,
    help='Path to bed file',
    required=False
)
@click.option(
    '--alpha',
    default = 2.0,
    help='alpha parameter',
    required=False
)
@click.option(
    '--beta',
    default = 1.0,
    help='beta parameter',
    required=False
)
@click.option(
    '--gamma',
    default = 1.0,
    help='gamma parameter',
    required=False
)
@click.option(
    '--pi_a',
    default = 0.3,
    help='pi_a parameter for HKY85 and F81 models',
    required=False
)
@click.option(
    '--pi_c',
    default = 0.2,
    help='pi_c parameter for HKY85 and F81 models',
    required=False
)
@click.option(
    '--pi_g',
    default = 0.2,
    help='pi_g parameter for HKY85 and F81 models',
    required=False
)
@click.option(
    '--pi_t',
    default = 0.3,
    help='pi_t parameter for HKY85 and F81 models',
    required=False
)
@click.option(
    '--omega',
    default = 0.5,
    help='omega parameter for K3P model',
    required=False
)
@click.option(
    '--model',
    '-m',
    type=click.Choice(['random', 'JC69', 'K2P', 'F81', 'HKY85', 'K3P']),
    default='JC69',
    help='evolutionary model to use for context calculation',
    required=True
)

# verbose flag
@click.option(
    '--verbose',
    '-V',
    is_flag=True,
    help='Print verbose output'
)
# flag to calculate dN/dS by gene
@click.option(
    '--by-read',
    '-r',
    is_flag=True,
    help='Calculate dN/dS by gene',
    required=False
)
def exhaust(fasta, bed, alpha, beta, gamma, pi_a, pi_c, pi_g, pi_t, omega, model, verbose, by_read):
    """
    Given a fasta file, calculate the dN/dS ratio using exhaustive method 
    where each permutation of the codon is tested
    """
    if by_read:
        if verbose:
            click.echo('Verbose mode enabled')
            click.echo('Exhaustive model started')
            click.echo(f"Calculating dN/dS ratio using {model} evolutionary model")
        dnds = dnds_exhaustive(fasta=fasta, bed=bed, alpha=alpha, beta=beta, gamma=gamma, pi_a=pi_a, pi_c=pi_c, pi_g=pi_g, pi_t=pi_t, omega=omega, model=model, by_read=True)
        if verbose:
            click.echo(f"dN/dS ratio for each gene: {dnds}")
        else:
            click.echo(f"{dnds}")
    else:
        if verbose:
            click.echo('Verbose mode enabled')
            click.echo('Exhaustive model started')
            click.echo(f"Calculating dN/dS ratio using {model} evolutionary model")
        dnds = dnds_exhaustive(fasta=fasta, bed=bed, alpha=alpha, beta=beta, gamma=gamma, pi_a=pi_a, pi_c=pi_c, pi_g=pi_g, pi_t=pi_t, omega=omega, model=model, by_read=False)
        if verbose:
            click.echo(f"dN/dS = {dnds}")
        else:
            click.echo(f"{dnds}")

if __name__ == '__main__':
    cli()
