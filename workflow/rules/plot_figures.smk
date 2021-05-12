## Plot the figures
rule plot_figures:
    input:
        "results/output/seasonalAbundance_BR.csv",
        "results/output/seasonalAbundance_NB.csv",
        "resources/species_names.csv",
        "results/output/empirical_connectivity_data.csv",
        "results/output/ORSIM_validation_metric.csv",
        "results/output/migratory_connectivity_strength.csv",
        "results/output/null_models_results.RData"
    output:
        "results/figures/Figure_2.png",
        "results/figures/Figure_3.pdf",
        "results/figures/Figure_S1.png",
        "results/figures/Figure_S5.png",
        "results/figures/Figure_S6.png",
        "results/figures/Figure_S7.png",
        "results/figures/Figure_S8.png",
        "results/figures/Figure_S9.png"
    conda:
        "../envs/env_R.yaml"
    script:
        "../scripts/20-figures/21-figures.R"
