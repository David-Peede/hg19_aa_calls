# Import packages.
import gzip

# Import custom library modules.
from hg19lib import load_hg19_seq

# Intialize the EPO directory.
EPO_DIR = '/oscar/data/ehuertas/dpeede/05_epo_aa_calls/hg19_aa_calls/data/homo_sapiens_ancestor_GRCh37_e71'



# Define a function to load a epo ancestral sequence.
def load_epo_seq(chrom):
    """
    Returns the EPO ancestral sequence for a focal chromosome.
    """
    # Intialize the ancestral sequence.
    anc_seq = ''
    # Open the fasta file.
    with open(f'{EPO_DIR}/homo_sapiens_ancestor_{chrom}.fa', 'r') as fa_data:
        # Compile the ancestral sequence.
        anc_seq = [line.strip().replace(' ', '') for line in fa_data if not line.startswith('>')]
    return ''.join(anc_seq)


# Define a function to build a EPO alignment dictionary.
def build_epo_aln_dicc(chrom):
    """
    Build the EPO alignment dictionary for a specified chromosome.
    """
    # Intialize a set of alleles to skip.
    skip_set = {'N', '-', '.'}
    # Intialize a dictionary to store the EPO alignment.
    epo_aln_dicc = {}
    # Load the Hg19 and EPO sequences.
    hg19_seq = load_hg19_seq(chrom)
    epo_seq = load_epo_seq(chrom)
    
    # Iterate through the chromosome.
    for pos, (hg19, epo) in enumerate(zip(hg19_seq, epo_seq), start=1):
        # If both alleles are valid.
        if hg19 not in skip_set and epo not in skip_set:
            # Update the dictionary.
            epo_aln_dicc[pos] = (hg19.upper(), epo.upper())
    return epo_aln_dicc