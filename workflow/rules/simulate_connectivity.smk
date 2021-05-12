## Download and process seasonal abundance estimates
rule seasonal_abundance:
    output:
        "results/output/seasonalAbundance_BR.csv",
        "results/output/seasonalAbundance_NB.csv",
        "results/output/distanceMatrix.csv"
    conda:
        "../envs/env_R.yaml"
    script:
        "../scripts/00-simulate_connectivity/01-seasonal_abundance.R"

## Run the optimal redistribution simulator
rule run_ORSIM:
    input:
        "results/output/seasonalAbundance_BR.csv",
        "results/output/seasonalAbundance_NB.csv",
        "results/output/distanceMatrix.csv"
    output:
        "results/output/ORSIM_results_woothr.csv"
    conda:
        "../envs/env_python.yaml"
    script:
        "../scripts/00-simulate_connectivity/02-orsim.py"
