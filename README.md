This repository contains the code to run the Optimal Redistribution Simulator that simulates animal migratory connectivity.

The code supports the analysis in Somveille M, Bay RA, Smith TB, Marra PP & Ruegg KC (2021) A general theory of avian migratory connectivity. Ecology Letters

The entire analysis can be reproduced using snakemake (see https://snakemake.readthedocs.io/en/stable/index.html for how to install and use snakemake). To run the analyis, the user can clone this repository onto their personal computer, open optimal-redistribution-simulator.Rproj in RStudio and in the terminal tab in RStudio run: snakemake --core 1. This will run the anlysis and generate the figures shown in Somveille et al. (2021) Ecology Letters, which are outputed in /results/figures.

Data is assumed to be stored in /resources.
