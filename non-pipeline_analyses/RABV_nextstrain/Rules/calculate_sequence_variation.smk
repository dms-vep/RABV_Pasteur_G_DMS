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
        protein_alignment = config["Ungapped_protein_alignment"],
    output:
        config["Glycoprotein_variation"],
    script:
        "../Scripts/calculate_variation.py"
