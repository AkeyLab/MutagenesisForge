import os
# return ratio of non-synonymous to synonymous substitutions

# give a measure of strength of natural selection
directory_path = ""

vep_files = os.listdir(directory_path)

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

