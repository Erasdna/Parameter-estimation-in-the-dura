#!/bin/bash -l

# Download the Gonzo dataset
wget "https://zenodo.org/records/14266867/files/mri-processed.zip?download=1" -O data/Gonzo.zip

# Unzip dataset
unzip data/Gonzo.zip -d data/Gonzo
# Remove large files we will not be using
rm data/Gonzo.zip
rm -r data/Gonzo/mri_dataset

#Run pipes:
uv run snakemake -s code/dura.smk --cores 1
uv run snakemake -s code/parameter_estimation.smk --cores "$1" --config qmc_samples="$2"
