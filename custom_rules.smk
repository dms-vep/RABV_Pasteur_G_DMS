"""Custom rules used in the ``snakemake`` pipeline.

This file is included by the pipeline ``Snakefile``.

"""


rule get_filtered_CSVs:
    """
    Get filtered DMS data CSVs.
    """
    input:
        func_effects="results/func_effects/averages/HEK293T_entry_func_effects.csv",
        func_effects_config="data/func_effects_config.yml",
        antibody_escape_config="data/antibody_escape_config.yml",
        nb="notebooks/get_filtered_CSVs.ipynb",
    params:
        antibody_escape_dir="results/antibody_escape/averages/",
        filtered_csv_dir="results/filtered_CSVs/",
    output:
        filtered_func_effects="results/filtered_CSVs/HEK293T_filitered_entry_func_effects.csv",
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
            -p antibody_escape_dir {params.antibody_escape_dir} \
            -p filtered_csv_dir {params.filtered_csv_dir} \
            -p filtered_func_effects {output.filtered_func_effects} \
            &> {log}
        """

rule compare_func_effects_and_natural_variation:
    """
    Compare DMS measured functional effects vs
    natural sequence variation among sequenced RABV G.
    """
    input:
        func_effects="results/filtered_CSVs/HEK293T_filitered_entry_func_effects.csv",
        func_effects_config="data/func_effects_config.yml",
        site_numbering_map="data/site_numbering_map.csv",
        natural_variation="non-pipeline_analyses/RABV_nextstrain/Results/G_variation.csv",
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
            -p site_numbering_map {input.site_numbering_map} \
            -p natural_variation {input.natural_variation} \
            -p natural_diversity_outdir {params.natural_diversity_outdir} \
            -p func_effects_vs_nature_html {output.func_effects_vs_nature_html} \
            &> {log}
        """

docs["Additional analyses and data files"] = {
    "Functional effects compared to natural sequence variation" : {
        "Functional effects vs natural diversity" : rules.compare_func_effects_and_natural_variation.output.func_effects_vs_nature_html,
        "Notebook comparing functional effects and natural variation" : rules.compare_func_effects_and_natural_variation.output.nb,
    },
    "Filtered CSVs for DMS data" : {
        "Notebook creating filtered CSVs for DMS data" : rules.get_filtered_CSVs.output.nb,
    },
}