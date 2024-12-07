# Description:
# Processes sequences using Augur

# Author:
# Caleb Carr

rule tree_sequences:
    """
    This rule creates a tree from sequences
    """
    input:
        alignment = "Results/{gene}/Alignments/codon_ungapped.fasta",
    output:
        tree = "Results/{gene}/Trees/tree_raw.nwk",
    params:
        outgroup = config["Outgroup"]
    conda:
        "../environment.yml",
    log:
        "Results/Logs/{gene}/tree_sequences.txt",
    shell:
        "augur tree "
        "--alignment {input.alignment} "
        "--tree-builder-args '\-o {params.outgroup}' "
        "--output {output.tree} &> {log}"


rule root_and_remove_outgroup:
    """
    This rule roots the tree and remove the outgroup
    """
    input:
        tree = "Results/{gene}/Trees/tree_raw.nwk",
        protein_alignment = "Results/{gene}/Alignments/protein_ungapped.fasta",
        codon_alignment = "Results/{gene}/Alignments/codon_ungapped.fasta",
    output:
        tree = "Results/{gene}/Trees/tree_no_outgroup.nwk",
        protein_alignment = "Results/{gene}/Alignments/protein_ungapped_no_outgroup.fasta",
        codon_alignment = "Results/{gene}/Alignments/codon_ungapped_no_outgroup.fasta",
    params:
        outgroup = config["Outgroup"],
    conda:
        "../environment.yml",
    log:
        "Results/Logs/{gene}/root_and_remove_outgroup.txt",
    script:
        "../Scripts/root_and_remove_outgroup.py"


rule refine_tree_sequences:
    """
    This rule refines the tree
    """
    input:
        alignment = "Results/{gene}/Alignments/codon_ungapped_no_outgroup.fasta",
        tree = "Results/{gene}/Trees/tree_no_outgroup.nwk",
        metadata = "Results/{gene}/metadata.tsv",
    output:
        tree = "Results/{gene}/Trees/tree.nwk",
        tree_nodes = "Results/{gene}/Trees/tree_branch_lengths.json",
    conda:
        "../environment.yml",
    log:
        "Results/Logs/{gene}/refine_tree_sequences.txt",
    shell:
        "augur refine "
        "--tree {input.tree} "
        "--alignment {input.alignment} "
        "--metadata {input.metadata} "
        "--output-tree {output.tree} "
        "--output-node-data {output.tree_nodes} "
        "--keep-root "
        "--verbosity 2 &> {log}"

rule traits_tree_sequences:
    """
    This rule gets traits from the tree
    """
    input:
        tree = "Results/{gene}/Trees/tree.nwk",
        metadata = "Results/{gene}/metadata.tsv",
    output:
        tree_traits = "Results/{gene}/Trees/tree_traits.json",
    conda:
        "../environment.yml",
    log:
        "Results/Logs/{gene}/traits_tree_sequences.txt",
    shell:
        "augur traits "
        "--tree {input.tree} "
        "--metadata {input.metadata} "
        "--output-node-data {output.tree_traits} " 
        "--columns region "
        "--confidence &> {log}"    
        
rule ancestral_tree_sequences:
    """
    This rule gets ancestral sequences from the tree
    """
    input:
        tree = "Results/{gene}/Trees/tree.nwk",
        alignment = "Results/{gene}/Alignments/codon_ungapped_no_outgroup.fasta",
    output:
        tree_muts = "Results/{gene}/Trees/tree_muts.json",
    conda:
        "../environment.yml",
    log:
        "Results/Logs/{gene}/ancestral_tree_sequences.txt",
    shell:
        "augur ancestral "
        "--tree {input.tree} "
        "--alignment {input.alignment} "
        "--output-node-data {output.tree_muts} " 
        "--inference joint &> {log}"

rule translate_tree_sequences:
    """
    This rule identifies amino-acid mutations from tree
    """
    input:
        tree = "Results/{gene}/Trees/tree.nwk",
        tree_muts = "Results/{gene}/Trees/tree_muts.json",
        reference = "Results/{gene}/{gene}_reference.gb",
    output:
        aa_muts = "Results/{gene}/Trees/tree_aa_muts.json",
        gene_alignments = "Results/{gene}/Augur_AA_Alignements/{gene}.fasta",
    conda:
        "../environment.yml",
    log:
        "Results/Logs/{gene}/translate_tree_sequences.txt",
    shell:
        "augur translate "
        "--tree {input.tree} "
        "--ancestral-sequences {input.tree_muts} "
        "--reference-sequence {input.reference} " 
        "--output-node-data {output.aa_muts} "
        "--genes {wildcards.gene} "
        "--alignment-output Results/{wildcards.gene}/Augur_AA_Alignements/%GENE.fasta &> {log}"


rule colors:
    """This rule assigns colors to traits"""
    input:
        color_schemes = config["Color_schemes"],
        color_orderings = config["Color_orderings"],
        metadata = "Results/{gene}/metadata.tsv",
    output:
        colors = "Results/{gene}/Trees/colors.tsv",
    conda:
        "../environment.yml",
    shell:
        "python Scripts/assign_colors.py "
        "--color-schemes {input.color_schemes} "
        "--ordering {input.color_orderings} "
        "--metadata {input.metadata} "
        "--output {output.colors}"


rule variant_escape_prediction:
    """This rule calculates variant escape"""
    input:
        escape_data = config["Antibody_escape_file"],
        alignment = "Results/{gene}/Alignments/protein_ungapped_no_outgroup.fasta",
        metadata = "Results/{gene}/metadata.tsv",
    params:
        antibody_name = lambda wildcards: config["Antibody_names"][wildcards.antibody],
        site_map = config["Site_map"],
        filter_params = config["Antibody_escape_prediction_params"]["Antibody_escape_filter_params"],
        reference_strain = config["Antibody_escape_prediction_params"]["Reference_DMS_strain"],
        site_column = config["Antibody_escape_prediction_params"]["Mutation_col"],
        escape_column = config["Antibody_escape_prediction_params"]["Mut_effect_col"],
    output:
        output_json = "Results/{gene}/Antibody_escape_predictions/{antibody}.json",
    conda:
        "../environment.yml",
    log:
        "Results/Logs/{gene}/{antibody}_variant_escape_prediction.txt",
    script:
        "../Scripts/assign_DMS_escape.py"
        

antibody_list = [
    "17C7",
    "CR4098",
    "CR57",
    "CTB012",
    "RVA122",
    "RVC20",
    "RVC58",
    "RVC68",
]

rule export_tree:
    """
    This rule exports the tree
    """
    input:
        tree = "Results/{gene}/Trees/tree.nwk",
        metadata = "Results/{gene}/metadata.tsv",
        tree_nodes = "Results/{gene}/Trees/tree_branch_lengths.json",
        tree_traits = "Results/{gene}/Trees/tree_traits.json",
        tree_muts = "Results/{gene}/Trees/tree_muts.json",
        aa_muts = "Results/{gene}/Trees/tree_aa_muts.json",
        auspice_config = config["Auspice_config"],
        colors = "Results/{gene}/Trees/colors.tsv",
        lat_longs = config["Lat_longs"],
        dms_predictions = expand("Results/{gene}/Antibody_escape_predictions/{antibody}.json", gene=["G"], antibody=antibody_list),
    params:
        title = lambda wildcards: config["Auspice_tree_titles"][wildcards.gene],
        description = lambda wildcards: config["Auspice_tree_descriptions"][wildcards.gene],
    output:
        auspice_tree = "auspice/{gene}.json",
    conda:
        "../environment.yml",
    shell:
        "augur export v2 "
        "--tree {input.tree} "
        "--title {params.title} "
        "--description {params.description} "
        "--metadata {input.metadata} "
        "--node-data {input.tree_nodes} {input.tree_traits} {input.tree_muts} {input.aa_muts} {input.dms_predictions} "
        "--include-root-sequence-inline "
        "--colors {input.colors} "
        "--lat-longs {input.lat_longs} "
        "--auspice-config {input.auspice_config} "
        "--output {output.auspice_tree}"