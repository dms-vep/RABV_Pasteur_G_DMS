"""Custom rules used in the ``snakemake`` pipeline.

This file is included by the pipeline ``Snakefile``.

"""


rule compare_func_effects_and_natural_variation:
    """
    Compare DMS measured functional effects vs
    natural sequence variation among sequenced RABV G.
    """
    input:
        func_effects="results/func_effects/averages/HEK293T_entry_func_effects.csv",
        site_numbering_map="data/site_numbering_map.csv",
        natural_variation="non-pipeline_analyses/RABV_nextstrain/Results/RABV_G_variation.csv",
        nb="notebooks/func_effects_vs_natural_variation.ipynb",
    params:
        min_times_seen=3,
    output:
        nb="results/notebooks/func_effects_vs_natural_variation.ipynb",
    conda:
        os.path.join(config["pipeline_path"], "environment.yml")
    log:
        "results/logs/compare_func_effects_and_natural_variation.txt",
    shell:
        """
        papermill {input.nb} {output.nb} \
            -p func_effects {input.func_effects} \
            -p site_numbering_map {input.site_numbering_map} \
            -p natural_variation {input.natural_variation} \
            -p min_times_seen {params.min_times_seen} \
            &> {log}
        """

docs["Additional analyses and data files"] = {
    "Functional effects compared to natural sequence variation" : {
        "Notebook comparing functional effects and natural variation" : rules.compare_func_effects_and_natural_variation.output.nb,
    }
}