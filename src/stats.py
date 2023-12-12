import os
import argparse
import glob

# return ratio of non-synonymous to synonymous substitutions

# give a measure of strength of natural selection
directory_path = os.getcwd()

pattern = 'vep*'

vep_files = glob.glob(os.path.join(directory_path, pattern))

def dnds(veps):
    syn = 0
    non_syn = 0
    for vep in veps:
        with open(vep, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                data = line.split()
                conseq = data[6]
                if conseq == 'synonymous_variant':
                    syn += 1
                else:
                    non_syn += 1
    ratio = non_syn / syn
    return(ratio)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument(
            '--stat', '-stat', type=str, required=True,
            help='statistical test to perform')
    
    args = parser.parse_args()

    stat_test = args.stat

    if stat_test == 'dnds':
        dnds(vep_files)

