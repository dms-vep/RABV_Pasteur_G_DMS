# Description:
# Aligns fasta sequences

# Author:
# Caleb Carr

rule align_protein_sequences:
    """
    This rule aligns all protein sequences using MAFFT
    """
    input:
        fasta_sequences = "Results/{gene}/ORF_sequences/protein.fasta",
    output:
        alignment = "Results/{gene}/Alignments/protein.fasta",
    conda:
        "../environment.yml",
    log:
        "Results/Logs/{gene}/align_protein_sequences.txt",
    shell:
        # --auto        Automatically selects an appropriate strategy from L-INS-i, 
        #               FFT-NS-i and FFT-NS-2, according to data size. Default: off (always FFT-NS-2)
        "mafft --auto {input.fasta_sequences} > {output.alignment} 2> {log}"

rule create_codon_alignment:
    """
    This rule creates a codon aligmnet from the
    codon sequences and protein alignment.
    """
    input: 
        protein_alignment = "Results/{gene}/Alignments/protein.fasta",
        codon_sequences = "Results/{gene}/ORF_sequences/codon.fasta",
    output:
        "Results/{gene}/Alignments/codon.fasta",
    conda:
        "../environment.yml",
    log:
        "Results/Logs/{gene}/create_codon_alignment.txt",
    script:
        "../Scripts/create_codon_alignment.py"


rule strip_protein_alignment_gaps:
    """
    This rule removes all gaps relative 
    to the reference sequence.
    """
    input: 
        alignment = "Results/{gene}/Alignments/protein.fasta",
    output:
        "Results/{gene}/Alignments/protein_ungapped.fasta",
    params:
        reference = config["Reference_accession"],
    conda:
        "../environment.yml",
    log:
        "Results/Logs/{gene}/strip_protein_alignment_gaps.txt",
    script:
        "../Scripts/remove_gaps_from_alignment.py"


rule strip_codon_alignment_gaps:
    """
    This rule removes all gaps relative 
    to the reference sequence.
    """
    input: 
        alignment = "Results/{gene}/Alignments/codon.fasta",
    output:
        "Results/{gene}/Alignments/codon_ungapped.fasta",
    params:
        reference = config["Reference_accession"],
    conda:
        "../environment.yml",
    log:
        "Results/Logs/{gene}/strip_codon_alignment_gaps.txt",
    script:
        "../Scripts/remove_gaps_from_alignment.py"