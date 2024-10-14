import pandas as pd
from io import StringIO

"""
this module contains functions for parsing and analyzing vep output
"""

def synoymous_detected(line):
    """Return True if the variant is synonymous."""
    df = pd.read_csv(StringIO(line), sep="\t")
    if df["Consequence"] == "synonymous_variant":
        return True
    return False


def total_dNds(veps):
    """Calculate dN/dS for a list of veps."""
    dN = 0
    dS = 0
    for vep in veps:
        for line in vep:
            if line.startswith("#"):
                continue
            if synoymous_detected(line):
                dS += 1
            else:
                dN += 1
    if dS == 0:
        return 0
    return dN / dS


def dNds(vep):
    """Calculate dN/dS for a single vep."""
    dN = 0
    dS = 0
    for line in vep:
        if line.startswith("#"):
            continue
        if synoymous_detected(line):
            dS += 1
        else:
            dN += 1
    if dS == 0:
        return 0
    return dN / dS
