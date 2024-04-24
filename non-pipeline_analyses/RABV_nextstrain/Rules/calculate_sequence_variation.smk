# Description:
# Processes and calculates sequence
# variation based on a MSA

# Author:
# Caleb Carr


rule remove_node_sequences:
    """
    This rule removes the tree node sequences
    from the MSA produced from Augur. 
    """
    input:
        raw_alignment = config["Gene_AA_alignments"],
    output:
        alignment = config["Glycoprotein_alignment"],
    conda:
        "../environment.yml",
    script:
        "../Scripts/process_alignment.py"


rule calculate_variation:
    """
    This rule calculates site level entropy and n effective 
    amino acids based on natural protein sequences.
    """
    input:
        protein_alignment = config["Glycoprotein_alignment"],
    output:
        config["Glycoprotein_variation"],
    script:
        "../Scripts/calculate_variation.py"
