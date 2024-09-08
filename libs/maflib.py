# Import packages.
import gzip



# Define a function to extract the sequences in an alignment block from a MAF file.
def maf_aln_block_generator(maf_file):
    """
    Iterates through a MAF file and yields the sequence alignment blocks.
    """
    # Intialize a list to store the sequences in the alignment block.
    aln_block = []
    # Iterate through every line.
    for line in maf_file:
        # If this is the start of a new alignment block.
        if line.startswith('a'):
            # If there are sequences from the previous alignment block.
            if aln_block:
                # Yield the sequences in the current alignment block.
                yield aln_block
                # Reset the list to store the sequences in the alignment block.
                aln_block = []
        # Else-if this a sequence within the alignment block.
        elif line.startswith('s'):
            # Update the alignment block.
            aln_block.append(line.strip())
    # If the last alignment block hasn't been processed.
    if aln_block:
        # Yield the sequences in the final alignment block.
        yield aln_block


# Define a function to compile the sequence fields for sequences in a MAF alignment block.
def compile_seq_field_info(maf_aln_block):
    """
    Extracts the sequence fields (see below) per sequence from a MAF alignment block.
    
    src: Name of the source sequence (alignment.chromosome).
    start: Start of the aligning region in the source sequence (zero-based).
    size: The length of the aligning region in the source sequence (non-dash characters only).
    strand: If "-", then the alignment is to the reverse-complemented source.
    srcSize: The size of the entire source sequence, not just the parts involved in the alignment.
    seq: The nucleotides in the alignment including any insertions ("-").
    """
    # Intialize a dictionary to store the sequence fields.
    seq_fields = {}
    # Iterate through every sequence in the alignment block.
    for seq_line in maf_aln_block:
        # Extract the sequence fields.
        src, start, size, strand, srcSize, seq = seq_line.split()[1:]
        # Update the sequence fields dictionary.
        seq_fields[src] = {
            'start': int(start), 'size': int(size),
            'strand': strand, 'srcSize': int(srcSize), 'seq': seq,
        }
    return seq_fields


# Define a function to build a MAF alignment dictionary.
def build_maf_aln_dicc(maf_file, chrom, aln_block_len_thresh=1):
    """
    Build a MAF alignment dictionary for a specified chromosome.
    
    maf_file: Bgziped .maf file.
    chrom: Chromosome to extract allele calls for.
    aln_block_len_thresh: Minimum alignment block length required to extract allele calls from.
    """
    # Intialize the chromosome id.
    chrom_id = f'chr{chrom}'
    # Intialize the reference sequence source.
    ref_seq_src = f'hg19.{chrom_id}'
    # Translation table to convert bases.
    vcf_trans = str.maketrans('-acgtn', 'NACGTN')
    # Intialize a dictionary to store the MAF alignment.
    maf_aln_dicc = {}
    
    # Open the MAF file.
    with gzip.open(maf_file, 'rt') as infile:
        # Iterate through every alignment block.
        for aln_block in maf_aln_block_generator(infile):
            # Extract the sequence fields.
            seq_fields = compile_seq_field_info(aln_block)
            # Grab the sequence fields for the refernce source.
            ref_seq_fields = seq_fields.get(ref_seq_src)
            
            # If the first sequence is the reference sequence and passes the alignment block length threshold.
            if (ref_seq_src == list(seq_fields.keys())[0]) and (seq_fields[ref_seq_src]['size'] >= aln_block_len_thresh):
                # Grab the sequence fields for the non-refernce source.
                alt_seq_fields = next((val for key, val in seq_fields.items() if key != ref_seq_src), None)
                
                # If the non-refernce source alignment is found.
                if alt_seq_fields:
                    # Grab the reference sequence.
                    ref_seq = ref_seq_fields['seq']
                    # Grab the length of the reference sequence (ie the number of ungapped bases).
                    ref_seq_len = ref_seq_fields['size']
                    # Intialize a generator for the non-gap indicies in the reference sequence alignment.
                    nogap_idx_iter = iter(i for i, base in enumerate(ref_seq) if base != '-')
                    # Generate the genomic positions.
                    positions = range(ref_seq_fields['start'] + 1, ref_seq_fields['start'] + 1 + ref_seq_len)
                    # Extract the the translated sequences.
                    ref_seq_trans = ref_seq.replace('-', '').translate(vcf_trans)
                    alt_seq_trans = alt_seq_fields['seq'].translate(vcf_trans)
                    
                    # Iterate through the reference sequence.
                    for j, pos in enumerate(positions):
                        # Update the dictionary.
                        maf_aln_dicc[pos] = (ref_seq_trans[j], alt_seq_trans[next(nogap_idx_iter)])
    return maf_aln_dicc