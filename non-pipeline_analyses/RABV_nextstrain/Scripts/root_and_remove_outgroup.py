# Description:
# Python script that roots a 
# phylogenetic tree on an outgroup
# and then removes that sequence

# Adapted from Jesse Bloom's root_and_remove_outgroup.py.ipynb

# Author:
# Caleb Carr

# Imports
import os
from Bio import SeqIO, Phylo

def root_and_remove_outgroup(input_tree_file, output_tree_file, outgroup):
    """Function to root tree and remove outgroup"""

    # Initialize tree
    tree = Phylo.read(input_tree_file, "newick")

    # Check that outgroup is in tree
    n_init = len(tree.get_terminals())
    assert any(clade.name == outgroup for clade in tree.get_terminals())

    # Root tree
    tree.root_with_outgroup(outgroup)
    tree.root = tree.root.clades[0]

    # Remove outgroup and check removal
    n_final = len(tree.get_terminals())
    assert not any(clade.name == outgroup for clade in tree.get_terminals())
    assert n_final == n_init - 1

    # Write new tree
    _ = Phylo.write(tree, output_tree_file, "newick")

def remove_sequence_from_alignment(input_alignment, output_alignment, outgroup):
    """Function to remove a sequence from an alignment"""

    # Initialize alignment and check sequence is in alignment
    alignment = list(SeqIO.parse(input_alignment, "fasta"))
    n_seqs_init = len(alignment)
    assert any(s.id == outgroup for s in alignment)

    # Remove sequence
    alignment = [s for s in alignment if s.id != outgroup]
    assert n_seqs_init == len(alignment) + 1

    # Write new alignment
    _ = SeqIO.write(alignment, output_alignment, "fasta")

def main():
    """
    Main method
    """

    # Input files
    input_tree = str(snakemake.input.tree)
    input_protein_alignment = str(snakemake.input.protein_alignment)
    input_codon_alignment = str(snakemake.input.codon_alignment)


    # Output files
    output_tree = str(snakemake.output.tree)
    output_protein_alignment = str(snakemake.output.protein_alignment)
    output_codon_alignment = str(snakemake.output.codon_alignment)
    output_log_file = str(snakemake.log)

    # Params
    outgroup = str(snakemake.params.outgroup)

    # Initialize log file
    log_file = open(output_log_file, "w")

    # Root tree without group and remove outgroup
    root_and_remove_outgroup(input_tree, output_tree, outgroup)

    # Remove outgroup from alignments
    remove_sequence_from_alignment(input_protein_alignment, output_protein_alignment, outgroup)
    remove_sequence_from_alignment(input_codon_alignment, output_codon_alignment, outgroup)

    log_file.write(f"Successfully rooted tree with {outgroup} and removed {outgroup} from alignments\n")

    # Close files
    log_file.close()


if __name__ == "__main__":
    main()