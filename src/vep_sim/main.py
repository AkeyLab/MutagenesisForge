from . import sim
from . import stat
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


@click.command()
# need to account for optional bed input
@click.argument(
    'input_bed_file',
    )
@click.argument(
    'input_mut_file', 
    )
@click.argument(
    'fasta_file', 
    )
@click.argument(
    'sim_num', 
    type = int, 
    )
@click.argument(
    'out',
    )

def main(input_bed_file, input_mut_file, fasta_file, sim_num, out):
    sim.sim(input_bed_file, input_mut_file, fasta_file, sim_num, out)

    # empirical vep output with yaml file
    emp_output = output_dir + "/" + "emp" + out + ".out"
    emp_cmd = emp_vep_path + emp_output
    emp_run = os.system(emp_cmd)

    # group vep output files
    sim_veps = find_files_with_prefix("sim")

    # dn/ds for grouped vep output files
    sim_total_dNds = stat.total_dNds(sim_veps)

    # empirical vep output files
    emp_veps = find_files_with_prefix("emp")
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