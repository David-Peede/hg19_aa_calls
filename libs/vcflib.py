# Import packages.
import gzip



# Define a function to identify what duplicated entires from a vcf.
def identify_dups(vcf, chrom):
    """
    Returns a set of duplicated positions.
    """
    # Intilailize the previous  position.
    prev_pos = None
    # Inialize a set to store duplicated positions.
    dup_set = set()
    # Read the vcf file.
    with gzip.open(vcf, 'rt') as data:
        # Iterate through every line in the original vcf file.
        for line in data:
            # If the line does not contains meta information.
            if not line.startswith('#'):
                # Split the line by tabs.
                spline = line.split()
                # If the line starts with the desired chromosome.
                if chrom == spline[0]:
                    # Grab the current poistion.
                    pos = spline[1]
                    # If the current position is already known to be duplicated.
                    if pos in dup_set:
                        # Update the previous position.
                        prev_pos = pos
                    # Else-if the current position is the same as the previous position.
                    elif pos == prev_pos:
                        # Add the position to the duplicated set.
                        dup_set.add(pos)
                    # Else, the current position is unique.
                    else:
                        # Update the previous position.
                        prev_pos = pos
    return dup_set