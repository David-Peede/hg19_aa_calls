# Import packages.
import argparse
from datetime import date
import gzip
import os
import sys

# Add the libraries directory to the system path for imports.
sys.path.append(os.path.abspath('/oscar/data/ehuertas/dpeede/05_epo_aa_calls/hg19_aa_calls/libs'))
# Import custom library modules.
from epolib import build_epo_aln_dicc



# Define a function extract the ancestral allele calls and output a VCF file.
def create_hg19_epo_vcf(chrom, buffer_size):
    """
    Writes a VCF file to standard out with the EPO ancestral allele for a specified chromosome.  
    """
    # Determine today's date.
    current_date = date.today().strftime('%m_%d_%Y')
    # Intialize the VCF header.
    vcf_header = f"""##fileformat=VCFv4.1
##FILTER=<ID=PASS,Description="All filters passed">
##source=create_hg19_epo_vcf_{current_date}
##contig=<ID={chrom}>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tEPO\n"""
    # Intialize a list to store the VCF lines.
    vcf_lines = [vcf_header]
    # Intialize a set of alleles to skip.
    skip_set = {'N', '-', '.'}
    # Build the EPO alignment dictionary.
    epo_aln = build_epo_aln_dicc(chrom)
    # Intialize the genotype conditions.
    is_chrom_y = chrom == 'Y'
    gt_same = '0' if is_chrom_y else '0|0'
    gt_diff = '1' if is_chrom_y else '1|1'
    
    # For every position.
    for pos in epo_aln:
        # Unpack the alleles.
        ref, alt = epo_aln[pos]
        # If both the reference and alternative allele are valid alleles.
        if ref not in skip_set and alt not in skip_set:
            # Determine the genotype.
            gt = gt_same if ref == alt else gt_diff
            # Update the VCF lines.
            vcf_lines.append(f'{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t.\tGT\t{gt}\n')
            # If the number of lines exceeds the buffer size.
            if len(vcf_lines) >= buffer_size:
                # Write the lines to the VCF file.
                sys.stdout.writelines(vcf_lines)
                # Clear the written lines.
                vcf_lines.clear()
    # If there are still lines to be written.
    if vcf_lines:
        # Write the remaining lines to the VCF file.
        sys.stdout.writelines(vcf_lines)
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


# Create the ancestral allele call VCF file.
create_hg19_epo_vcf(args.chrom, args.buffer_size)