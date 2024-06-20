# Description:
# Processes sequences using Augur

# Author:
# Caleb Carr

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
        "../dmsa-pred/dmsa_env.yaml"
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
        "../dmsa-pred/dmsa_env.yaml"
    shell:
        "augur traits "
        "--tree {input.tree} "
        "--metadata {input.metadata} "
        "--output-node-data {output.tree_traits} " 
        "--columns country phylogroup region "
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
        "../dmsa-pred/dmsa_env.yaml"
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
        "../dmsa-pred/dmsa_env.yaml"
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
    conda:
        "../dmsa-pred/dmsa_env.yaml"
    shell:
        "python Scripts/assign_colors.py "
        "--color-schemes {input.color_schemes} "
        "--ordering {input.color_orderings} "
        "--metadata {input.metadata} "
        "--output {output.colors}"


rule variant_escape_prediction:
    """This rule calculates variant escape"""
    input:
        alignment = config["Gene_AA_alignments"],
    output:
        node_data = "Results/dmsa-phenotype/{collection}/{experiment}_escape_prediction.json",
        pred_data = "Results/dmsa-phenotype/{collection}/{experiment}_escape_prediction.csv"
    params:
        basedir = lambda w: config["dmsa_phenotype_collections"].get(w.collection)['mut_effects_dir'],
        dms_wt_seq_id = "NC_001542_2018-08-13",
        mut_effect_col = lambda w: config["dmsa_phenotype_collections"].get(w.collection)['mut_effect_col'],
        mutation_col = lambda w: config["dmsa_phenotype_collections"].get(w.collection)['mutation_col'],
        mut_effects_df = lambda w: os.path.join(
            config["dmsa_phenotype_collections"].get(w.collection)['mut_effects_dir'], 
            w.experiment
        ),
    conda:
        "../dmsa-pred/dmsa_env.yaml"
    shell:
        "python dmsa-pred/dmsa_pred.py phenotype-prediction "
        "--model-type additive "
        "--alignment {input.alignment} "
        "--dms-wt-seq-id {params.dms_wt_seq_id} "
        "--mask-seqs-with-disallowed-aa-subs False "
        "--min-pred-pheno 0.0 "
        "--mut-effects-df {params.mut_effects_df} "
        "--mut-effect-col {params.mut_effect_col} "
        "--mutation-col {params.mutation_col} "
        "--experiment-label {wildcards.experiment} "
        "--output-json {output.node_data} "
        "--output-df {output.pred_data}"


def _get_variant_escape_node_data(wildcards):
    inputs=[]
    wildcards_dict = dict(wildcards)

    import glob
    for collection_name, collection_dict in config['dmsa_phenotype_collections'].items():

        # run the predictions using every csv in the glob
        requested_files = expand(
            rules.variant_escape_prediction.output.node_data,
            collection=collection_name,
            experiment=[
                os.path.basename(fp) 
                for fp in glob.glob(collection_dict['mut_effects_dir']+"/*.csv")
            ],
            **wildcards_dict
        )
        inputs.extend(requested_files)

    return inputs


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
        lat_longs = config["Lat_longs"],
        escape_predictions = _get_variant_escape_node_data,
    output:
        auspice_tree = config["Auspice_tree"],
    conda:
        "../dmsa-pred/dmsa_env.yaml"
    shell:
        "augur export v2 "
        "--tree {input.tree} "
        "--metadata {input.metadata} "
        "--node-data {input.tree_nodes} {input.tree_traits} {input.tree_muts} {input.aa_muts} {input.escape_predictions} "
        "--include-root-sequence "
        "--colors {input.colors} "
        "--lat-longs {input.lat_longs} "
        "--auspice-config {input.auspice_config} "
        "--output {output.auspice_tree}"