# scdna_replication_paper

Code for performing analysis associated with Weiner et al., "Inferring replication timing and proliferation dynamics from single-cell DNA sequencing data", Nature Communications (2024). 

Snakemake pipelines are used to process data, run the [PERT model for scRT inference](https://github.com/shahcompbio/scdna_replication_tools), and do compute-intensive analysis tasks. For each unique dataset (i.e. `fitness`, `sig_lines`, etc), there is a snakemake rule file at `rules/*.smk` and a set of analysis scripts in the `scripts/*/` folder. The output from these snakemake pipelines are used to create all main figures `notebooks/main_figures/` and supplementary figures `notebooks/supplementary_figures/`. For a detailed instructions on how to use PERT for your own data, see the [tutorials in the scdna_replication_tools github repository](https://github.com/shahcompbio/scdna_replication_tools/tree/main/notebooks).

All source data needed to reproduce the figures (inputs to the jupyter notebooks), along with detailed descriptions of the file structure, are provided in our corresponding [Zenodo submission](https://doi.org/10.5281/zenodo.12786373). Please note that you will need to change input and output file paths if you wish to use this code and online source data to reproduce any of our analysis.
