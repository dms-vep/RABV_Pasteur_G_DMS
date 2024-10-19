# Description:
# Processes nucleotide fastas to find largest ORF in 
# forward direction and translates to protein sequences

# Author:
# Caleb Carr


rule get_protein_sequences:
    """
    This rule runs EMBOSS getorf to extract 
    amino acid sequences. 
    """
    input:
        sequences = "Results/{gene}/nucleotide.fasta",
    params:
        min_ORF_size = lambda wildcards: config["ORF_min_size"][wildcards.gene],
    output:
        "Results/{gene}/ORF_sequences/protein_temp.fasta",
    conda:
        "../environment.yml",
    log:
        "Results/Logs/{gene}/get_protein_sequences.txt",
    shell:
        # The '-sequence' flag signals for the input file
        # while the '-outseq' file signals for the output 
        # file. The '-find 1' flag means the amino acid 
        # sequences between START and STOP codons are returned.
        # The '-minsize' flag means minimum nucleotide size of 
        # ORF to report. The '-reverse' flag signals if ORFs on 
        # the reverse strand should be found as well. 
        "getorf -sequence {input.sequences} -outseq {output} -find 1 -minsize {params.min_ORF_size} -reverse No -verbose &> {log}"


rule get_codon_sequences:
    """
    This rule runs EMBOSS getorf to extract 
    codon sequences. 
    """
    input:
        sequences = "Results/{gene}/nucleotide.fasta",
    params:
        min_ORF_size = lambda wildcards: config["ORF_min_size"][wildcards.gene],
    output:
        "Results/{gene}/ORF_sequences/codon_temp.fasta",
    conda:
        "../environment.yml",
    log:
        "Results/Logs/{gene}/get_codon_sequences.txt",
    shell:
        # The '-sequence' flag signals for the input file
        # while the '-outseq' file signals for the output 
        # file. The '-find 3' flag means the nucleotide 
        # sequences between START and STOP codons are returned.
        # The '-minsize' flag means minimum nucleotide size of 
        # ORF to report. The '-reverse' flag signals if ORFs on 
        # the reverse strand should be found as well. 
        "getorf -sequence {input.sequences} -outseq {output} -find 3 -minsize {params.min_ORF_size} -reverse No -verbose &> {log}"


rule check_number_ORFs_found:
    """
    This rule checks both the amino acid
    and codon fastas to verify there is a 
    one-to-one mapping from nucleotide
    seqeunce to protein and codon sequences. 
    """
    input:
        nucleotide = "Results/{gene}/nucleotide.fasta",
        protein = "Results/{gene}/ORF_sequences/protein_temp.fasta",
        codon = "Results/{gene}/ORF_sequences/codon_temp.fasta",
    output:
        "Results/{gene}/ORF_sequences/ORF_extraction_verification.txt"
    conda:
        "../environment.yml",
    script:
        "../Scripts/extracted_sequence_verification.py"


rule edit_protein_fasta_headers:
    """
    This rule edits the protein sequence fasta
    headers to remove the appended info from EMBOSS
    getorf.
    """
    input: 
        sequences = "Results/{gene}/ORF_sequences/protein_temp.fasta",
    output:
        sequences = "Results/{gene}/ORF_sequences/protein.fasta",
    conda:
        "../environment.yml",
    log:
        "Results/Logs/{gene}/edit_protein_fasta_headers.txt",
    script:
        "../Scripts/edit_fasta_headers.py"


rule edit_codon_fasta_headers:
    """
    This rule edits the protein sequence fasta
    headers to remove the appended info from EMBOSS
    getorf.
    """
    input: 
        sequences = "Results/{gene}/ORF_sequences/codon_temp.fasta",
    output:
        sequences = "Results/{gene}/ORF_sequences/codon.fasta",
    conda:
        "../environment.yml",
    log:
        "Results/Logs/{gene}/edit_codon_fasta_headers.txt",
    script:
        "../Scripts/edit_fasta_headers.py"