# Description:
# Processes and calculates sequence
# variation based on a MSA

# Author:
# Caleb Carr


rule calculate_variation:
    """
    This rule calculates site level entropy and n effective 
    amino acids based on natural protein sequences.
    """
    input:
        protein_alignment = "Results/{gene}/Alignments/protein_ungapped_no_outgroup.fasta",
    output:
        "Results/{gene}/Alignments/{gene}_natural_variation.csv",
    conda:
        "../environment.yml",
    script:
        "../Scripts/calculate_variation.py"
