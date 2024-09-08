# Import packages.
import gzip

# Intialize the Hg19 directory.
HG19_DIR = '/oscar/data/ehuertas/dpeede/05_epo_aa_calls/hg19_aa_calls/data/hg19'


# Define a function the load the hg19 reference sequence.
def load_hg19_seq(chrom):
    """
    Returns the Hg19 reference sequence for a focal chromosome.
    """
    # Intialize a flag for compiling the reference sequence.
    compile_seq = False
    # Intialize a list for storing reference sequence.
    ref_seq = []
    # Open the fasta file.
    with gzip.open(f'{HG19_DIR}/hg19.fa.gz', 'rt') as fa_data:
        # For every line in the fasta file.
        for line in fa_data:
            # Clean the line.
            line = line.strip()
            # If this is a header line.
            if line.startswith('>'):
                # Check to see if this is the focal chromosome.
                if line.lstrip('>') == f'chr{chrom}':
                    # Start compiling the reference sequence.
                    compile_seq = True
                # Else-if this is the next chromosome's header line. 
                elif compile_seq:
                    # Stop iterating after compiling the refrence sequence.
                    break
            # Else-if we are iterating through the focal chromosome.
            elif compile_seq:
                    # Update the reference sequence.
                    ref_seq.append(line)
    return ''.join(ref_seq)