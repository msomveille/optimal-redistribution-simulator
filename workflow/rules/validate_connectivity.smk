## Prepare data for model validation
rule data_prep:
    input:
        "results/output/seasonalAbundance_BR.csv",
        "results/output/seasonalAbundance_NB.csv",
        "resources/species_names.csv"
    output:
        "results/output/empirical_connectivity_data.csv"
    conda:
        "../envs/env_R.yaml"
    script:
        "../scripts/10-validate_connectivity/11-data_prep.R"

## Compute validation metric for each species
rule validation_metric:
    input:
        "results/output/seasonalAbundance_BR.csv",
        "results/output/seasonalAbundance_NB.csv",
        "resources/species_names.csv",
        "results/output/empirical_connectivity_data.csv"
    output:
        "results/output/ORSIM_validation_metric.csv"
    conda:
        "../envs/env_R.yaml"
    script:
        "../scripts/10-validate_connectivity/12-validation_metric.R"

## Run null models
rule null_models:
    input:
        "results/output/seasonalAbundance_BR.csv",
        "results/output/seasonalAbundance_NB.csv",
        "resources/species_names.csv",
        "results/output/empirical_connectivity_data.csv",
        "results/output/ORSIM_validation_metric.csv"
    output:
        "results/output/null_models_results.RData"
    conda:
        "../envs/env_R.yaml"
    script:
        "../scripts/10-validate_connectivity/13-null_models.R"

## Calculate empirical and simulated strenght of migratory conectivity
rule conectivity_strength:
    input:
        "results/output/seasonalAbundance_BR.csv",
        "results/output/seasonalAbundance_NB.csv",
        "resources/species_names.csv",
        "results/output/empirical_connectivity_data.csv"
    output:
        "results/output/migratory_connectivity_strength.csv"
    conda:
        "../envs/env_R.yaml"
    script:
        "../scripts/10-validate_connectivity/14-conectivity_strength.R"
