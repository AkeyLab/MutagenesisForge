"""
this module contain functions for calculating transition/transversion ratio
"""

def tstv(vcf):

    transition = ["AG", "GA", "CT", "TC"]
    transversion = ["AC", "CA", "AT", "TA", "GC", "CG", "GT", "TG"]

    transition_count = 0
    transversion_count = 0

    with open(vcf, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            line = line.strip().split()
            ref_base = line[3]
            alt = line[4]
            mut = ref_base + alt
            if mut in transition:
                transition_count += 1
            if mut in transversion:
                transversion_count += 1
            else:
                continue
    return transition_count / transversion_count
