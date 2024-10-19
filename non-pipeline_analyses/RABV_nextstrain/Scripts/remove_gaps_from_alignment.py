# Description:
# Python script that removes gaps 
# from an alignment relative to a 
# reference sequence

# Author:
# Caleb Carr

# Imports
from Bio import AlignIO

# Functions
def remove_gaps_from_alignment(alignment_file, reference, ungapped_alignment, output_log_file):

    # Initialize alignment and reference
    alignment = AlignIO.read(alignment_file, "fasta")
    reference_sequence = None

    # Initialize log file
    log_file = open(output_log_file, "w")
    
    # Get reference sequence
    for seq in alignment:
        if reference in seq.description:
            reference_sequence = str(seq.seq)

    # Check that reference was found
    assert reference_sequence != None, "ERROR: Reference sequence is not in alignment!"

    # Keep track of gaps already moved
    num_gaps_removed = 0
    for index,char in enumerate(reference_sequence):
        # Remove if gap found 
        if char == "-":
            alignment = alignment[:, :(index-num_gaps_removed)] + alignment[:, (index-num_gaps_removed+1):]
            num_gaps_removed += 1


    # Write new alignment file
    AlignIO.write(alignment, ungapped_alignment, "fasta")

    log_file.write("Success!")

    # Close files
    log_file.close()





def main():
    """
    Main method
    """

    # Input files
    alignment_file = str(snakemake.input.alignment)

    # Params
    reference = str(snakemake.params.reference)

    # Output files
    ungapped_alignment = str(snakemake.output)
    output_log_file = str(snakemake.log)

    remove_gaps_from_alignment(alignment_file, reference, ungapped_alignment, output_log_file)


if __name__ == "__main__":
    main()