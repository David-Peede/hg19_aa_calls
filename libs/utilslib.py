# Define a function to find the all, shared, and private positions between two alignment dictionaries.
def compile_position_sets(aln_dicc1, aln_dicc2):
    """
    Find the union, intersection, and set differences of positions between two alignment dictionaries.
    """
    # Extract the positions.
    pos1 = set(aln_dicc1)
    pos2 = set(aln_dicc2)
    # Determine and sort the union.
    all_pos = sorted(pos1 | pos2)
    # Determine the intersection.
    shared_pos = pos1 & pos2
    # Determine the set differences.
    priv_pos1 = pos1 - pos2
    priv_pos2 = pos2 - pos1
    return all_pos, shared_pos, priv_pos1, priv_pos2

# Define a function to build the meta information dictionary for the tgp.
def build_tgp_meta_dicc(tgp_meta_file):
    """
    Returns a dictionary where the keys are an individual's index in the VCF
    file and the values are an individual's superpopulation and population IDs.
    """
    # Intialize a dictionary.
    ind_dicc = {}
    # Open the meta info file.
    with open(tgp_meta_file, 'r') as pop_data:
        # Skip the header line.
        next(pop_data)
        # Intialize a column counter.
        c_col = 9
        # For every line.
        for i, line in enumerate(pop_data):
            # Split the line.
            spline = line.split()
            # Grab the current population and superpopulation.
            pop, spop = spline[1], spline[2]
            # Update the dictionary.
            ind_dicc[c_col + i] = (spop, pop)
    return ind_dicc