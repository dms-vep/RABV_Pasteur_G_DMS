# Description:
# Python script to process the 
# raw gene alignments produced 
# from Augur

# Author:
# Caleb Carr

# Imports
import datetime
import os
import pandas as pd
from Bio import Entrez, SeqIO, AlignIO
Entrez.email = "RABVexample@RABVexample.com"     # Always tell NCBI who you are

# Functions
def remove_node_sequences(raw_alignment, alignment):

    total_seq_count = 0
    retained_seq_count = 0
    # Iterate through alignments and remove node sequenes
    with open(alignment, "w") as new_protein_alignment:
        for curr_fasta in AlignIO.read(raw_alignment, "fasta"):
            total_seq_count += 1
            curr_name = str(curr_fasta.description)
            if "NODE_" in curr_name:
                continue
            retained_seq_count += 1
            SeqIO.write(curr_fasta, new_protein_alignment, "fasta")

    print(f"{total_seq_count} total sequences processed and {retained_seq_count} sequences retained")

    # Close files
    new_protein_alignment.close()



def main():
    """
    Main method
    """

    # Input files
    raw_alignment = str(snakemake.input.raw_alignment) 

    # Output files
    alignment = str(snakemake.output.alignment) 

    remove_node_sequences(raw_alignment, alignment)


if __name__ == "__main__":
    main()