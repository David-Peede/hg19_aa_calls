# Import packages.
import argparse
import gzip
import os
import sys

# Add the libraries directory to the system path for imports.
sys.path.append(os.path.abspath('/oscar/data/ehuertas/dpeede/05_epo_aa_calls/hg19_aa_calls/libs'))
# Import custom library modules.
from utilslib import build_tgp_meta_dicc
from vcflib import identify_dups

# Intialize the directory paths.
DATA_DIR = '/oscar/data/ehuertas/data/modern_genomes/1000genomes/ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502'
TABLE_DIR = '/oscar/data/ehuertas/dpeede/05_epo_aa_calls/hg19_aa_calls/tgp_info'

# Intialize global variables.
TGP_SPOP_LIST = ['AFR', 'SAS', 'EAS', 'EUR', 'AMR']
TGP_POP_LIST = [
    'LWK', 'GWD', 'MSL', 'ESN', 'YRI', 'ACB', 'ASW', # AFR.
    'BEB', 'STU', 'ITU', 'PJL', 'GIH', # SAS.
    'CHB', 'KHV', 'CHS', 'JPT', 'CDX', # EAS.    
    'TSI', 'CEU', 'IBS', 'GBR', 'FIN', # EUR.
    'PEL', 'MXL', 'CLM', 'PUR', # AMR.
]
TGP_META_FILE = f'{DATA_DIR}/integrated_call_samples_v3.20130502.ALL.panel'



# Define a function to qc a line from the unfiltered all archaics merged vcf.
def process_tgp_autosome_line(line, chrom, dup_set, ind_dicc):
    """
    Returns the qc information for an autosomal TGP VCF line.
    """
    # Split the line by tabs.
    spline = line.split()
    # Grab the position, rsID, refernce, and alternative alleles.
    pos, rs_id, ref, alt = spline[1:5]
    
    # If the current position has duplicated records.
    if pos in dup_set:
        # Return the qc information.
        return f'{chrom}\t{pos}\t{rs_id}\t{ref}\t{alt}\tdup_record\n', None
    # Else-if the length of refernce and alternative fields is not exactly 2,
    # ie if the site is a not mono- or biallelic snp.
    elif (len(ref) + len(alt)) != 2:
        # Return the qc information.
        return f'{chrom}\t{pos}\t{rs_id}\t{ref}\t{alt}\tmulti_allelic_or_sv\n', None

    # Else, the position is a bi-allelic snp.
    else:
        # Intialize superpopulation and population dictionaries.
        spop_dicc = {spop: {'aac': 0, 'chr': 0} for spop in TGP_SPOP_LIST}
        pop_dicc = {pop: {'aac': 0, 'chr': 0} for pop in TGP_POP_LIST}
        # For every individual.
        for i, (spop, pop) in ind_dicc.items():
            # Count the number of alternative alleles.
            aac = spline[i].count('1')
            # Update the superpopulation and population counts.
            spop_dicc[spop]['aac'] += aac
            spop_dicc[spop]['chr'] += 2
            pop_dicc[pop]['aac'] += aac
            pop_dicc[pop]['chr'] += 2
        # Compile the alternative allele frequencies.
        spop_aaf = [spop_dicc[spop]['aac'] / float(spop_dicc[spop]['chr']) for spop in TGP_SPOP_LIST]
        pop_aaf = [pop_dicc[pop]['aac'] / float(pop_dicc[pop]['chr']) for pop in TGP_POP_LIST]
        spop_pop_aaf_str = '\t'.join(map(str, spop_aaf + pop_aaf))
        # Return the qc information.
        return None, f'{chrom}\t{pos}\t{rs_id}\t{ref}\t{alt}\t{spop_pop_aaf_str}\n'

# Define a function to parse an autosomal tgp vcf file and output a qc and alternative allele frequency table.
def create_tgp_autosomal_tables(vcf, chrom, tgp_meta_file, buffer_size):
    """
    Creates two gzipped TXT files:
    
    1) tgp_failed_qc_report_chr{chrom}.txt.gz - Information for the sites that failed QC.
    2) tgp_spop_pop_aaf_biallelic_snps_chr{chrom}.txt.gz - Alternative allele frequencies
       for every superpopulation and population at bi-allelic SNPs.
    """
    # Load the meta information dictionary.
    ind_dicc = build_tgp_meta_dicc(tgp_meta_file)
    # Determine the duplicated positions.
    dup_set = identify_dups(vcf, chrom)
    # Intialize a list to store the the qc information.
    failed_qc_lines = ['CHROM\tPOS\trsID\tREF\tALT\tFAILED_CONDITION\n']
    spop_pop_cols = '\t'.join(TGP_SPOP_LIST + TGP_POP_LIST)
    passed_qc_lines = [f'CHROM\tPOS\trsID\tREF\tALT\t{spop_pop_cols}\n']
    
    # Read the vcf file and intilialize the qc reports.
    with gzip.open(vcf, 'rt') as vcf_data, \
         gzip.open(f'{TABLE_DIR}/qc_tables/tgp_failed_qc_report_chr{chrom}.txt.gz', 'wt') as failed_qc_file, \
         gzip.open(f'{TABLE_DIR}/aaf_tables/tgp_spop_pop_aaf_biallelic_snps_chr{chrom}.txt.gz', 'wt') as passed_qc_file:
        # Iterate through every line in the original vcf file.
        for line in vcf_data:
            # If the line does not contains meta information.
            if not line.startswith('#'):
                # Process the current line.
                failed_qc_info, passed_qc_info = process_tgp_autosome_line(line, chrom, dup_set, ind_dicc)
                # If the current line failed qc.
                if failed_qc_info:
                    # Append the qc file list.
                    failed_qc_lines.append(failed_qc_info)
                # Else, the current line passed qc.
                else:
                    # Append the qc file list.
                    passed_qc_lines.append(passed_qc_info)
            # If the number of failed qc lines exceeds the buffer size.
            if len(failed_qc_lines) >= buffer_size:
                # Write the failed qc lines to the file.
                failed_qc_file.writelines(failed_qc_lines)
                # Clear the written qc lines.
                failed_qc_lines.clear()
            # If the number of passed qc lines exceeds the buffer size.
            if len(passed_qc_lines) >= buffer_size:
                # Write the passed qc lines to the file.
                passed_qc_file.writelines(passed_qc_lines)
                # Clear the written qc lines.
                passed_qc_lines.clear()

        # If there are still failed qc lines to be written.
        if failed_qc_lines:
            # Write the remaining failed qc lines.
            failed_qc_file.writelines(failed_qc_lines)
        # If there are still passed qc lines to be written.
        if passed_qc_lines:
            # Write the remaining passed qc lines.
            passed_qc_file.writelines(passed_qc_lines)
    return


# Intialize the command line options.
parser = argparse.ArgumentParser()
parser.add_argument(
    '-c', '--chrom', required=True, type=str, action='store',
    help='Chromosome ID of the TGP VCF file.',
)
parser.add_argument(
    '-b', '--buffer_size', required=False, type=int, action='store', default=100_000,
    help='Maximum number of lines to be stored in the buffer before writting (default = 100,000 lines).',
)
args = parser.parse_args()

# Intialize the vcf file.
VCF_FILE = f'{DATA_DIR}/ALL.chr{args.chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz'
# Create the qc and alterantive allele frequency tables.
create_tgp_autosomal_tables(VCF_FILE, args.chrom, TGP_META_FILE, args.buffer_size)