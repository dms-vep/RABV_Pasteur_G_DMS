"""Custom rules used in the ``snakemake`` pipeline.

This file is included by the pipeline ``Snakefile``.

"""


with open(config["dms_viz_config"]) as f:
    dms_viz_config = yaml.safe_load(f)

rule configure_dms_viz:
    """Configure a JSON file for `dms-viz`."""
    input:
        unpack(
            lambda wc: (
                {"per_antibody_escape_csv": 
                    "results/summaries/all_antibodies_and_cell_entry_per_antibody_escape.csv"
                }
                if dms_viz_config[wc.struct]["escape_or_phenotype"] == "escape"
                else {}
            )
        ),
        phenotypes_csv="results/summaries/all_antibodies_and_cell_entry.csv",
    output:
        json="results/dms-viz/{struct}.json",
        json_no_description=temp("results/dms-viz/{struct}_no_description.json"),
        sitemap=temp("results/dms-viz/{struct}_sitemap.csv"),
        phenotypes=temp("results/dms-viz/{struct}_phenotypes.csv"),
        description_md=temp("results/dms-viz/{struct}_description.md"),
        pdb_file=temp("results/dms-viz/{struct}.pdb"),
    params:
        config=lambda wc: dms_viz_config[wc.struct],
        description_suffix=(
            f"Study by {config['authors']} ({config['year']}).\n"
            f"See [{config['github_repo_url']}]({config['github_repo_url']}) for code."
        ),
    log:
        notebook="results/notebooks/configure_dms_viz_{struct}.ipynb",
    conda:
        "envs/dms-viz.yml"
    notebook:
        "notebooks/configure_dms_viz.py.ipynb"

docs["Visualizations of DMS data on protein structure"] = {
    "dms-viz JSONs": {
        "JSON files for dms-viz": {
            struct: rules.configure_dms_viz.output.json.format(struct=struct)
            for struct in dms_viz_config
        },
        "Notebooks generating JSONs": {
            struct: rules.configure_dms_viz.log.notebook.format(struct=struct)
            for struct in dms_viz_config
        },
    },
}


rule get_filtered_CSVs:
    """
    Get filtered DMS data CSVs.
    """
    input:
        # dummy files to confirm its run in correct order
        mAb_RVC20="results/antibody_escape/averages/RVC20_mut_effect.csv",
        mAb_RVA122="results/antibody_escape/averages/RVA122_mut_effect.csv",
        mAb_17C7="results/antibody_escape/averages/17C7_mut_effect.csv",
        mAb_RVC58="results/antibody_escape/averages/RVC58_mut_effect.csv",
        mAb_CR4098="results/antibody_escape/averages/CR4098_mut_effect.csv",
        mAb_CR57="results/antibody_escape/averages/CR57_mut_effect.csv",
        mAb_CTB012="results/antibody_escape/averages/CTB012_mut_effect.csv",
        mAb_RVC68="results/antibody_escape/averages/RVC68_mut_effect.csv",
        func_effects="results/func_effects/averages/HEK293T_entry_func_effects.csv",
        func_effects_config="data/func_effects_config.yml",
        antibody_escape_config="data/antibody_escape_config.yml",
        site_numbering_map="data/site_numbering_map.csv",
        nb="notebooks/get_filtered_CSVs.ipynb",
    params:
        antibody_escape_dir="results/antibody_escape/averages/",
        filtered_antibody_csv_dir="results/filtered_antibody_escape_CSVs/",
        filtered_func_csv_dir="results/filtered_func_effects_CSV/",
    output:
        filtered_func_effects="results/filtered_func_effects_CSV/HEK293T_filitered_entry_func_effects.csv",
        nb="results/notebooks/get_filtered_CSVs.ipynb",
    conda:
        os.path.join(config["pipeline_path"], "environment.yml")
    log:
        "results/logs/get_filtered_CSVs.txt",
    shell:
        """
        papermill {input.nb} {output.nb} \
            -p func_effects {input.func_effects} \
            -p func_effects_config {input.func_effects_config} \
            -p antibody_escape_config {input.antibody_escape_config} \
            -p site_numbering_map {input.site_numbering_map} \
            -p antibody_escape_dir {params.antibody_escape_dir} \
            -p filtered_antibody_csv_dir {params.filtered_antibody_csv_dir} \
            -p filtered_func_csv_dir {params.filtered_func_csv_dir} \
            -p filtered_func_effects {output.filtered_func_effects} \
            &> {log}
        """

rule compare_func_effects_and_natural_variation:
    """
    Compare DMS measured functional effects vs
    natural sequence variation among sequenced RABV G.
    """
    input:
        func_effects="results/filtered_func_effects_CSV/HEK293T_filitered_entry_func_effects.csv",
        func_effects_config="data/func_effects_config.yml",
        natural_variation="non-pipeline_analyses/RABV_nextstrain/Results/G/Alignments/G_natural_variation.csv",
        nb="notebooks/func_effects_vs_natural_variation.ipynb",
    params:
        natural_diversity_outdir="results/natural_diversity_comparison/",
    output:
        func_effects_vs_nature_html="results/natural_diversity_comparison/func_effects_vs_natural_diversity.html",
        nb="results/notebooks/func_effects_vs_natural_variation.ipynb",
    conda:
        os.path.join(config["pipeline_path"], "environment.yml")
    log:
        "results/logs/compare_func_effects_and_natural_variation.txt",
    shell:
        """
        papermill {input.nb} {output.nb} \
            -p func_effects {input.func_effects} \
            -p func_effects_config {input.func_effects_config} \
            -p natural_variation {input.natural_variation} \
            -p natural_diversity_outdir {params.natural_diversity_outdir} \
            -p func_effects_vs_nature_html {output.func_effects_vs_nature_html} \
            &> {log}
        """


rule compare_antibody_escape_and_natural_variation:
    """
    Compare DMS measured antibody escape vs
    natural sequence variation among sequenced RABV G.
    """
    input:
        # dummy file to confirm its run in correct order
        filtered_func_effects="results/filtered_func_effects_CSV/HEK293T_filitered_entry_func_effects.csv",
        sequence_metadata="non-pipeline_analyses/RABV_nextstrain/Results/G/metadata.tsv",
        sequence_alignment="non-pipeline_analyses/RABV_nextstrain/Results/G/Alignments/protein_ungapped_no_outgroup.fasta",
        nb="notebooks/natural_sequence_antibody_escape.ipynb",
    params:
        filtered_antibody_csv_dir="results/filtered_antibody_escape_CSVs/",
        natural_diversity_outdir="results/natural_diversity_comparison/",
    output:
        mutation_count_csv="results/natural_diversity_comparison/mutation_counts.csv",
        antibody_escape_vs_nature_html="results/natural_diversity_comparison/antibody_escape_vs_natural_diversity.html",
        nb="results/notebooks/natural_sequence_antibody_escape.ipynb",
    conda:
        os.path.join(config["pipeline_path"], "environment.yml")
    log:
        "results/logs/compare_antibody_escape_and_natural_variation.txt",
    shell:
        """
        papermill {input.nb} {output.nb} \
            -p sequence_metadata {input.sequence_metadata} \
            -p sequence_alignment {input.sequence_alignment} \
            -p filtered_antibody_csv_dir {params.filtered_antibody_csv_dir} \
            -p natural_diversity_outdir {params.natural_diversity_outdir} \
            -p mutation_count_csv {output.mutation_count_csv} \
            -p antibody_escape_vs_nature_html {output.antibody_escape_vs_nature_html} \
            &> {log}
        """


docs["Additional analyses and data files"] = {
    "Comparison to natural sequence variation" : {
        "Functional effects vs natural diversity" : rules.compare_func_effects_and_natural_variation.output.func_effects_vs_nature_html,
        "Antibody escape vs natural diversity" : rules.compare_antibody_escape_and_natural_variation.output.antibody_escape_vs_nature_html,
    },
    "Filtered CSVs for DMS data" : {
        "Notebook creating filtered CSVs for DMS data" : rules.get_filtered_CSVs.output.nb,
    },
}
