# Description:
# Aligns fasta sequences

# Author:
# Caleb Carr

rule align_protein_sequences:
    """
    This rule aligns all GPC protein sequences using MAFFT
    """
    input:
        fasta_sequences = config["Protein_sequences"],
    output:
        alignment = config["Protein_alignment"],
    shell:
        # --auto        Automatically selects an appropriate strategy from L-INS-i, 
        #               FFT-NS-i and FFT-NS-2, according to data size. Default: off (always FFT-NS-2)
        "mafft --auto {input.fasta_sequences} > {output.alignment}"

rule create_codon_alignment:
    """
    This rule creates a codon aligmnet from the
    codon sequences and protein alignment.
    """
    input: 
        protein_alignment = config["Protein_alignment"],
        codon_sequences = config["Codon_sequences"],
    output:
        config["Codon_alignment"]
    script:
        "../Scripts/create_codon_alignment.py"


rule strip_protein_alignment_gaps:
    """
    This rule removes all gaps relative 
    to the reference sequence.
    """
    input: 
        codon_alignment = config["Protein_alignment"],
    output:
        config["Ungapped_protein_alignment"]
    params:
        reference = config["Reference_accession"],
    script:
        "../Scripts/remove_gaps_from_alignment.py"


rule strip_codon_alignment_gaps:
    """
    This rule removes all gaps relative 
    to the reference sequence.
    """
    input: 
        codon_alignment = config["Codon_alignment"],
    output:
        config["Ungapped_codon_alignment"]
    params:
        reference = config["Reference_accession"],
    script:
        "../Scripts/remove_gaps_from_alignment.py"