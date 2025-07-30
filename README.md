# miniLyndon experiments


All the programs needed to reproduce the experiments are managed with [conda](https://github.com/conda/conda).
It was tested with version 25.5.1, but similar versions should work.
Compatible packages managers (e.g., mamba or micromamba) should work too.

The experiments were executed and tested under Ubuntu 24.04.

Commands needed to reproduce the experiments:

```sh
conda env create -p ./local-env -f environment.yml
conda activate ./local-env
snakemake --sdm=conda --conda-create-envs-only

## Execute all the experiments
snakemake --sdm=conda --cores=16   ## Replace 16 with the number of CPU cores you want to use

## Prepare the MultiQC reports
snakemake --cores=16 --sdm=conda multiqc_all

## Compute statistics about the simulated reads
snakemake --cores=16 --sdm=conda read_stats_all

```
