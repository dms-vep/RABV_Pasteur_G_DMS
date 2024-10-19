# Description:
# Python script that removes the extra
# info appended to fasta headers from EMBOSS
# getorf

# Author:
# Caleb Carr

# Imports
import os
from Bio import SeqIO

# Functions
def edit_fasta_headers(input_sequences, output_sequences, output_log_file):

    # Initialize log file
    log_file = open(output_log_file, "w")

    # Iterate through headers and remove added info
    with open(output_sequences, "w") as new_sequences:
        for curr_fasta in SeqIO.parse(input_sequences, "fasta"):
            new_name = str(curr_fasta.description).split(" ")[0][:-2]
            log_file.write(f"{str(curr_fasta.description)} processed to ")
            curr_fasta.id = new_name
            curr_fasta.description = ""
            SeqIO.write(curr_fasta, new_sequences, "fasta")
            log_file.write(f"{new_name}\n")

    # Close files
    new_sequences.close()
    log_file.close()

def main():
    """
    Main method
    """

    # Input files
    input_sequences = str(snakemake.input.sequences)

    # Output files
    output_sequences = str(snakemake.output.sequences)
    output_log_file = str(snakemake.log)

    edit_fasta_headers(input_sequences, output_sequences, output_log_file)


if __name__ == "__main__":
    main()