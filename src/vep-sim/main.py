from . import sim
from . import stat
import os
import click
import numpy as np
from scipy import stats
import yaml

# load parameters
def load_params():
    with open("params.yaml", "r") as f:
        params = yaml.safe_load(f)
    return params

parameters = load_params('parameters.yaml')

# access parameters
emp_vep_path = parameters.get('emp_vep_path')


@click.command()
@click.option('--input_bed_file', help='Input bed file')
@click.option('--input_mut_file', help='Input mut file')
@click.option('--fasta_file', help='Input fasta file')
@click.option('--sim_num', type = int, default = 10, help='Number of simulations')

def main(input_bed_file, input_mut_file, fasta_file, sim_num):
    sim.sim(input_bed_file, input_mut_file, fasta_file, sim_num)

    # empirical vep output with yaml file
    emp_run = os.system(vep_path)

    # group vep output files
    sim_veps = [] 

    # dn/ds for grouped vep output files
    sim_total_dNds = stat.total_dNds(sim_veps)

    # empirical vep output files
    emp_veps = []
    emp_total_dNds = stat.total_dNds(emp_veps)

    # individual vep output files dn/ds
    ind_sim_dNds = []
    for vep in sim_veps:
        ind_sim_dNds.append(stat.dNds(vep))

    # statistical analysis
    sim_std_dev = np.std(ind_sim_dNds)
    
    t_statistic, p_value = stats.ttest_1samp(emp_total_dNds, ind_sim_dNds)

    # write data to terminal
    print(f"Empirical dN/dS: {emp_total_dNds}")
    print(f"Simulated dN/dS: {sim_total_dNds}")
    print(f"Standard deviation: {sim_std_dev}")
    print(f"t-statistic: {t_statistic}")
    print(f"p-value: {p_value}")

    
if __name__ == '__main__':
    main()