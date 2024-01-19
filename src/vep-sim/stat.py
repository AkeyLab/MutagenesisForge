def total_dNds(veps):
    """Calculate dN/dS for a list of veps."""
    dN = 0
    dS = 0
    for vep in veps:
        if vep['Consequence'] == 'missense_variant':
            dN += 1
        elif vep['Consequence'] == 'synonymous_variant':
            dS += 1
    if dS == 0:
        return 0
    return dN/dS

def dNds(vep):
    """Calculate dN/dS for a single vep."""
    dN = 0
    dS = 0
    if vep['Consequence'] == 'missense_variant':
        dN += 1
    elif vep['Consequence'] == 'synonymous_variant':
        dS += 1
    if dS == 0:
        return 0
    return dN/dS
