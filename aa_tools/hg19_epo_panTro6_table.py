# Import packages.
import argparse
import gzip
import os
import sys

# Add the libraries directory to the system path for imports.
sys.path.append(os.path.abspath('/oscar/data/ehuertas/dpeede/05_epo_aa_calls/hg19_aa_calls/libs'))
# Import custom library modules.
from hg19lib import load_hg19_seq
from epolib import build_epo_aln_dicc
from maflib import build_maf_aln_dicc
from utilslib import compile_position_sets

# Intialize the directory paths.
DATA_DIR = '/oscar/data/ehuertas/dpeede/05_epo_aa_calls/hg19_aa_calls/data'
TABLE_DIR = '/oscar/data/ehuertas/dpeede/05_epo_aa_calls/hg19_aa_calls/aa_calls/tables'



# Define a function to construct a table to compare two alignments.
def create_epo_vs_panTro6_table(chrom, buffer_size):
    """
    Creates a gzipped CSV file comparing the EPO vs. panTro6 ancestral allele for a specified chromosome.  
    """
    # Intialize a list to store the table information.
    table_info = ['CHROM,POS,Hg19,EPO,panTro6\n']
    # Load the Hg19 sequence.
    hg19_seq = load_hg19_seq(chrom).upper()
    # Build the EPO and MAF alignment dictionaries.
    epo_aln = build_epo_aln_dicc(chrom)
    panTro6_aln = build_maf_aln_dicc(f'{DATA_DIR}/panTro/hg19.panTro6.synNet.maf.gz', chrom)
    # Compile the position information.
    union_pos, intersect_pos, uniq_epo_pos, uniq_panTro6_pos = compile_position_sets(epo_aln, panTro6_aln)
    
    # Open the table file.
    with gzip.open(f'{TABLE_DIR}/hg19_epo_panTro6_chr{chrom}.csv.gz', 'wt') as table_file:
        # Iterate through all the positions between the two alignments.
        for pos in union_pos:
            # If the position is in both alignments.
            if pos in intersect_pos:
                # Update the table information.
                table_info.append(f'{chrom},{pos},{hg19_seq[pos - 1]},{epo_aln[pos][1]},{panTro6_aln[pos][1]}\n')
            # Else-if the position is unique to the EPO alignment.
            elif pos in uniq_epo_pos:
                # Update the table information.
                table_info.append(f'{chrom},{pos},{hg19_seq[pos - 1]},{epo_aln[pos][1]},N\n')
            # Else-if the position is unique to the panTro6 alignment.
            elif pos in uniq_panTro6_pos:
                # Update the table information.
                table_info.append(f'{chrom},{pos},{hg19_seq[pos - 1]},N,{panTro6_aln[pos][1]}\n')
            
            # If the number of lines exceeds the buffer size.
            if len(table_info) >= buffer_size:
                # Write the lines to the table file.
                table_file.writelines(table_info)
                # Clear the written lines.
                table_info.clear()
        # If there are still lines to be written.
        if table_info:
            # Write the remaining lines to the table file.
            table_file.writelines(table_info)
    return


# Intialize the command line options.
parser = argparse.ArgumentParser()
parser.add_argument(
    '-c', '--chrom', required=True, type=str, action='store',
    help='Chromosome to construct the ancestral allele table for.',
)
parser.add_argument(
    '-b', '--buffer_size', required=False, type=int, action='store', default=100_000,
    help='Maximum number of lines to be stored in the buffer before writting (default = 100,000 lines).',
)
args = parser.parse_args()


# Create the ancestral allele call comparison table.
create_epo_vs_panTro6_table(args.chrom, args.buffer_size)