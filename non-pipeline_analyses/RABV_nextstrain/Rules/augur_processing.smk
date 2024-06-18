# Description:
# Processes sequences using Augur

# Author:
# Caleb Carr

# rule index_sequences:
#     """
#     This rule indexes the sequences
#     """
#     input:
#         sequences = config["Nucleotide_unfiltered_sequences"],
#     output:
#         sequence_index = config["Nucleotide_sequences_indexed"],
#     conda:
#         "../environment.yml",
#     shell:
#         "augur index "
#         "--sequences {input.sequences} "
#         "--output {output.sequence_index}"

# rule filter_sequences:
#     """
#     This rule filters sequences
#     """
#     input:
#         fasta_sequences = config["Nucleotide_unfiltered_sequences"],
#         sequence_index = config["Nucleotide_sequences_indexed"],
#         metadata = config["Metadata"],
#         excluded_strains = config["Dropped_strains"],
#     output:
#         fasta_sequences = config["Nucleotide_sequences"],
#     conda:
#         "../environment.yml",
#     shell:
#         "augur filter "
#         "--sequences {input.fasta_sequences} "
#         "--sequence-index {input.sequence_index} "
#         "--metadata {input.metadata} "
#         "--exclude {input.excluded_strains} "
#         "--output {output.fasta_sequences}"

# rule align_sequences:
#     """
#     This rule aligns the sequences
#     """
#     input:
#         sequences = config["Nucleotide_sequences"],
#         reference = config["Reference_genbank"],
#     output:
#         alignment = config["Nucleotide_alignment"],
#     conda:
#         "../environment.yml",
#     shell:
#         "augur align "
#         "--sequences {input.sequences} "
#         "--reference-sequence {input.reference} "
#         "--output {output.alignment} "
#         "--remove-reference "
#         "--fill-gaps"

rule tree_sequences:
    """
    This rule creates a tree from sequences
    """
    input:
        alignment = config["Ungapped_codon_alignment"],
    output:
        tree = config["Nucleotide_raw_tree"],
    conda:
        "../environment.yml",
    shell:
        "augur tree "
        "--alignment {input.alignment} "
        "--output {output.tree}"

rule refine_tree_sequences:
    """
    This rule refines the tree
    """
    input:
        alignment = config["Ungapped_codon_alignment"],
        tree = config["Nucleotide_raw_tree"],
        metadata = config["Metadata"],
    output:
        tree = config["Nucleotide_tree"],
        tree_nodes = config["Nucleotide_tree_nodes"]
    conda:
        "../environment.yml",
    shell:
        "augur refine "
        "--tree {input.tree} "
        "--alignment {input.alignment} "
        "--metadata {input.metadata} "
        "--output-tree {output.tree} "
        "--output-node-data {output.tree_nodes} "
        "--timetree "
        "--coalescent opt "
        "--date-confidence "
        "--date-inference marginal"
        # "--clock-filter-iqd 4"

rule traits_tree_sequences:
    """
    This rule gets traits from the tree
    """
    input:
        tree = config["Nucleotide_tree"],
        metadata = config["Metadata"],
    output:
        tree_traits = config["Nucleotide_tree_traits"]
    conda:
        "../environment.yml",
    shell:
        "augur traits "
        "--tree {input.tree} "
        "--metadata {input.metadata} "
        "--output-node-data {output.tree_traits} " 
        "--columns country phylogroup "
        "--confidence"    
        
rule ancestral_tree_sequences:
    """
    This rule gets ancestral sequences from the tree
    """
    input:
        tree = config["Nucleotide_tree"],
        alignment = config["Ungapped_codon_alignment"],
    output:
        tree_muts = config["Nucleotide_tree_mutations"]
    conda:
        "../environment.yml",
    shell:
        "augur ancestral "
        "--tree {input.tree} "
        "--alignment {input.alignment} "
        "--output-node-data {output.tree_muts} " 
        "--inference joint"

rule translate_tree_sequences:
    """
    This rule identifies amino-acid mutations from tree
    """
    input:
        tree = config["Nucleotide_tree"],
        tree_muts = config["Nucleotide_tree_mutations"],
        reference = config["Reference_genbank"],
    output:
        aa_muts = config["Nucleotide_AA_mutations"],
        gene_alignments = expand("Results/Augur_AA_Alignements/{gene}.fasta", gene=["Glycoprotein"]),
    conda:
        "../environment.yml",
    shell:
        "augur translate "
        "--tree {input.tree} "
        "--ancestral-sequences {input.tree_muts} "
        "--reference-sequence {input.reference} " 
        "--output-node-data {output.aa_muts} "
        "--genes Glycoprotein "
        "--alignment-output Results/Augur_AA_Alignements/%GENE.fasta"


rule colors:
    """This rule assigns colors to traits"""
    input:
        color_schemes = config["Color_schemes"],
        color_orderings = config["Color_orderings"],
        metadata = config["Metadata"],
    output:
        colors = config["Colors"],
    shell:
        "python Scripts/assign_colors.py "
        "--color-schemes {input.color_schemes} "
        "--ordering {input.color_orderings} "
        "--metadata {input.metadata} "
        "--output {output.colors}"


rule export_tree:
    """
    This rule exports the tree
    """
    input:
        tree = config["Nucleotide_tree"],
        metadata = config["Metadata"],
        tree_nodes = config["Nucleotide_tree_nodes"],
        tree_traits = config["Nucleotide_tree_traits"],
        tree_muts = config["Nucleotide_tree_mutations"],
        aa_muts = config["Nucleotide_AA_mutations"],
        auspice_config = config["Auspice_config"],
        colors = config["Colors"],
    output:
        auspice_tree = config["Auspice_tree"],
    conda:
        "../environment.yml",
    shell:
        "augur export v2 "
        "--tree {input.tree} "
        "--metadata {input.metadata} "
        "--node-data {input.tree_nodes} {input.tree_traits} {input.tree_muts} {input.aa_muts} "
        "--include-root-sequence "
        "--colors {input.colors} "
        "--auspice-config {input.auspice_config} "
        "--output {output.auspice_tree}"