import collections
# import gzip potentially not needed
import pysam
# import pandas as pd potentially not needed
import numpy as np

def exhaustive(path, by_read=False):
    # need to overhall this to support fasta files not in exact format of the one I downloaded from NCBI

    # codon to amino acid table
    codon_to_amino = {
        'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
        'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
        'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
        'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
        'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
        'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
        'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
        'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
    }

    # Precalculate the synonymous, missense, and nonsense counts for all mutations of each codon at each position
    synonymous_muts_per_codon = {}
    missense_muts_per_codon = {}
    nonsense_muts_per_codon = {}

    bases = {'A','C','G','T'}
    positions = [0, 1, 2]

    #iterate through all 64 codons
    for codon in codon_to_amino:
        num_synonymous_muts = 0
        num_missense_muts = 0
        num_nonsense_muts = 0
        
        #iterate through all 3 positions in each codon
        for pos in positions:
            base = codon[pos]
            amino = codon_to_amino[codon]
            possible_mutations = bases.difference(base)
            
            #iterate through all 3 possible mutations of the base at the current position
            for mutated_base in possible_mutations:
                mutated_codon = codon[:pos]+mutated_base+codon[pos+1:]
                mutated_amino = codon_to_amino[mutated_codon]
                
                #assign the mutation to be nonsense/missense/synonymous
                if mutated_amino == '*':
                    num_nonsense_muts += 1
                elif amino != mutated_amino:
                    num_missense_muts += 1
                else:
                    num_synonymous_muts += 1
                    
                    
                #print(f'orig_codon={codon}, mutated_codon={mutated_codon}, orig_amino={amino}, mut_amino={mutated_amino}')
        
        #store the number of each type of mutation for this codon
        synonymous_muts_per_codon[codon] = num_synonymous_muts
        missense_muts_per_codon[codon] = num_missense_muts
        nonsense_muts_per_codon[codon] = num_nonsense_muts
            
    fasta = pysam.FastaFile(path)

    data = {
        'has_ATG_start':[],
        'stop_codon_end':[],
        'cds_length':[],
        'num_codons':[],
        'num_synonymous':[],
        'num_missense':[],
        'num_nonsense':[],
    }

    for ref in fasta.references:
        # fetching gene coding region data
        seq = fasta.fetch(ref)
        #print(seq + '\n' + '\n')  #debugging
        
        # string manipulation to isolate end and start of cds
        #cds = [s for s in ref_data if s.startswith("CDS:")]
        #start, end = cds[0].split(':')[1].split('-')
        # type converstion to int
        #start = int(start)
        #end = int(end)
        #cds_seq = seq[start-1:end]
        
        #sanity checks
        #- the cDNA seq length should divisible by 3
        #- the first codon should be ATG
        #- the last should be a stop codon (TAA, TGA, TAG)
        #assert len(cds_seq)%3 == 0                 #actually it seems like this isn't always true somehow
        #assert cds_seq[:3] == 'ATG'                #same with this
        #assert cds_seq[-3:] in {'TAA','TGA','TAG'} #same with this
        
        has_atg_start = False
        if seq[:3] == "ATG":
            has_atg_start = True
        
        stop_codons = ["TAA", "TGA", "TAG"]
        stop_codon_end = False
        if seq[-3:] in stop_codons:
            stop_codon_end = True
        
        
        #count how many of each type of codon there are
        codon_frequency = collections.Counter(seq[i:i+3] for i in range(0, len(seq), 3))
        
        #remove "codons" that are less than 3 basepairs (I don't know why the transcript lens aren't all div by 3)
        codon_frequency = {codon:count for codon,count in codon_frequency.items() if len(codon) == 3}
        
        #use the pre-calculated counts of synonymous/missense/nonsense per codon to count total mutations
        num_synonymous = sum(
            synonymous_muts_per_codon[codon]*count
            for codon,count in codon_frequency.items()
        )
        
        num_missense = sum(
            missense_muts_per_codon[codon]*count
            for codon,count in codon_frequency.items()
        )

        num_nonsense = sum(
            nonsense_muts_per_codon[codon]*count
            for codon,count in codon_frequency.items()
        )

        #update the "data" dictionary
        #do not need to store ensg, enst, gene_name
        data['has_ATG_start'].append(has_atg_start)
        data['stop_codon_end'].append(stop_codon_end)
        data['num_codons'].append(sum(codon_frequency.values()))
        data['num_synonymous'].append(num_synonymous)
        data['num_missense'].append(num_missense)
        data['num_nonsense'].append(num_nonsense)


    #close the fasta file
    fasta.close()


    #convert the "data" dictionary into dnds
    #use the file-wide methodology
    dnds_method_1 = (sum(data['num_missense']) + sum(data['num_nonsense'])) / sum(data['num_synonymous'])
    dnds_method_2 = []
    for i in range(len(data['has_ATG_start'])):
        if data['num_synonymous'][i] == 0:
            continue
        dnds_method_2.append((data['num_missense'][i] + data['num_nonsense'][i]) / data['num_synonymous'][i])
    
    dnds2_mean = np.mean(dnds_method_2)

    #convert the "data" dictionary to be a pandas table
    #df = pd.DataFrame(data)
    #df['dnds'] = (df['num_nonsense'] + df['num_missense'])/df['num_synonymous']
    #dnds_values = df.iloc[-1].tolist()
    if by_read:
        return(dnds2_mean)
    
    return(dnds_method_1)


