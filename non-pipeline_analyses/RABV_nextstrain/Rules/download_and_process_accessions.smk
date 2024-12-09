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
        max_frac_N = config["max_frac_N"]["genes"],
        accesstions_to_exclude = config["Accessions_to_exclude"],
        accessions_to_include = config["Accessions_to_include"],
        remove_duplicates = config["Remove_duplicates"],
        asia = config["Asia"],
        oceania = config["Oceania"],
        africa = config["Africa"],
        europe = config["Europe"],
        south_america = config["South America"],
        north_america = config["North America"],
        genbank_features = config["Genbank_features"],
        host_grouper = config["Host_grouper"],
        phylogroups = config["Phylogroups"],
        gene_names = lambda wildcards: config["Gene_groupings"][wildcards.gene],
        desired_gene = lambda wildcards: wildcards.gene,
    output:
        fasta_sequences = "Results/{gene}/nucleotide.fasta",
        metadata = "Results/{gene}/metadata.tsv",
    conda:
        "../environment.yml",
    log:
        "Results/Logs/{gene}/download_and_process_accessions.txt",
    script:
        "../Scripts/download_NCBI_sequences.py"


rule create_new_references:
    """
    This rule runs a script to create new references.
    """
    input:
        config["Reference_genbank"],
    params:
        gene_name = lambda wildcards: config["Gene_groupings"][wildcards.gene][0],
    output:
        new_reference = "Results/{gene}/{gene}_reference.gb",
    conda:
        "../environment.yml",
    log:
        "Results/Logs/{gene}/create_new_references.txt",
    shell:
        "python Scripts/newreference.py "
        "--reference {input} "
        "--output-genbank {output.new_reference} "
        "--new-name G_full_length "
        "--gene {params.gene_name} &> {log}"
