# Description:
# Processes nucleotide fastas to find largest ORF in 
# forward direction and translates to protein sequences

# Author:
# Caleb Carr


rule get_protein_sequences:
    """
    This rule runs EMBOSS getorf to extract 
    amino acid sequences for GPC segment. 
    """
    input:
        sequences = config["Nucleotide_sequences"],
    params:
        config["Min_ORF_threshold"],
    output:
        config["Protein_sequences_temp"],
    conda:
        "../environment.yml",
    shell:
        # The '-sequence' flag signals for the input file
        # while the '-outseq' file signals for the output 
        # file. The '-find 1' flag means the amino acid 
        # sequences between START and STOP codons are returned.
        # The '-minsize' flag means minimum nucleotide size of 
        # ORF to report. The '-reverse' flag signals if ORFs on 
        # the reverse strand should be found as well. 
        "getorf -sequence {input.sequences} -outseq {output} -find 1 -minsize {params} -reverse No -verbose"


rule get_codon_sequences:
    """
    This rule runs EMBOSS getorf to extract 
    codon sequences for G. 
    """
    input:
        sequences = config["Nucleotide_sequences"],
    params:
        config["Min_ORF_threshold"],
    output:
        config["Codon_sequences_temp"],
    conda:
        "../environment.yml",
    shell:
        # The '-sequence' flag signals for the input file
        # while the '-outseq' file signals for the output 
        # file. The '-find 3' flag means the nucleotide 
        # sequences between START and STOP codons are returned.
        # The '-minsize' flag means minimum nucleotide size of 
        # ORF to report. The '-reverse' flag signals if ORFs on 
        # the reverse strand should be found as well. 
        "getorf -sequence {input.sequences} -outseq {output} -find 3 -minsize {params} -reverse No -verbose"


rule check_number_ORFs_found:
    """
    This rule checks both the amino acid
    and codon fastas to verify there is a 
    one-to-one mapping from G nucleotide
    seqeunce to protein and codon sequences. 
    """
    input:
        nucleotide = config["Nucleotide_sequences"],
        protein = config["Protein_sequences_temp"],
        codon = config["Codon_sequences_temp"],
    output:
        config["Sequence_verification_log"]
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
        sequences = config["Protein_sequences_temp"],
    output:
        sequences = config["Protein_sequences"],
    conda:
        "../environment.yml",
    script:
        "../Scripts/edit_fasta_headers.py"


rule edit_codon_fasta_headers:
    """
    This rule edits the protein sequence fasta
    headers to remove the appended info from EMBOSS
    getorf.
    """
    input: 
        sequences = config["Codon_sequences_temp"],
    output:
        sequences = config["Codon_sequences"],
    conda:
        "../environment.yml",
    script:
        "../Scripts/edit_fasta_headers.py"