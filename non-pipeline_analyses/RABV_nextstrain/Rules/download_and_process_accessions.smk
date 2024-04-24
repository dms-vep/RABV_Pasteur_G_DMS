# Description:
# Downloads accessions and process/extract sequences and metadata

# Author:
# Caleb Carr


rule download_and_process_accessions:
    """
    This rule runs the download_NCBI_sequences.py script which
    downloads all genbank files from a list of accessions and 
    extracts metadata. 
    """
    input:
        accession_list = config["Accession_list"],
    params:
        genome_size_threshold_lower = config["Genome_size_threshold_lower"],
        genome_size_threshold_upper = config["Genome_size_threshold_upper"],
        max_frac_N = config["max_frac_N"],
        desired_segment = "",
        accesstions_to_exclude = config["Accessions_to_exclude"],
    output:
        fasta_sequences = config["Nucleotide_unfiltered_sequences"],
        metadata = config["Metadata"],
    conda:
        "../environment.yml",
    script:
        "../Scripts/download_NCBI_sequences.py"
