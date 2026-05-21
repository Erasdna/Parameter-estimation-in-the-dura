# Multiphasic Transport in the Human Cranial Dura: How Meningeal Lymphatic Vessels Drive Waste Clearance from Cerebrospinal Fluid

This repository contains supplementary code for the paper
> Girelli et al. 2026.
> Multiphasic Transport in the Human Brain Dura: How Meningeal Lymphatic Vessels Drive Waste Clearance from Brain Cerebrospinal Fluid, arxiv


# Abstract

**Background:**
  -
  The discovery of meningeal lymphatic vessels (mLVs) has reshaped our understanding of brain waste clearance, positioning them as a very important clearance route. However, quantifying the precise hydrodynamics and solute transport mechanisms within the human dura mater \emph{in vivo }remains highly challenging due to the multiscale nature of cranial fluid dynamics and imaging resolution limits.

**Methods:**
  -
  We integrated patient-specific gMRI data with a novel anatomical segmentation of the cranial dura mater. By employing a hierarchy of inverse mathematical and micro-mechanical models, we estimated previously inaccessible physiological transport parameters. These baseline estimates subsequently informed a spatially distributed, double-porosity 3D computational model of the human dura to simulate the fluid dynamics of cranial clearance.

**Results:**
-
Our analysis suggests that purely diffusive macroscopic models fail to replicate in vivo tracer distribution. Conversely, a dual-porosity advection-diffusion framework successfully captures the tissue kinetics. The transport within the mLV compartment is characterized by highly convective clearance (with $Pe>1$), with macroscopic fluid velocities and effective lymphatic diffusion coefficients confirming an advective capacity sufficient to clear more than the brain's total daily protein production.

**Conclusion:**
-
Solute clearance within the human cranial dura mater is a multiphasic, advection-dominated process. Rather than relying on passive Fickian diffusion, molecules are actively and rapidly cleared through the meningeal lymphatic network via convective mechanisms. This hierarchical modeling approach bridges clinical imaging and microscopic fluid dynamics, clarifying the primary mechanical drivers of brain homeostasis.


## Reproducing

### Getting started

First install `uv`
```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```
Then create a virtual environment and install dependencies

```bash
uv venv
source .venv/bin/activate
uv sync --locked --all-extras --dev
```



The [Gonzo dataset](https://zenodo.org/records/14266867) can be downloaded directly from Zenodo:
```bash
  wget "https://zenodo.org/records/14266867/files/mri-processed.zip?download=1" -O data/Gonzo.zip
  unzip data/Gonzo.zip -d data/Gonzo
```

Or alternatively using [mri-toolkit](https://github.com/scientificcomputing/mri-toolkit)
```bash
  uv run mritk datasets download gonzo -o data/Gonzo
  mv data/Gonzo/mri-processed/mri_processed_data data/Gonzo/mri_processed_data
```
This second option will also download various files from the dataset which are not used in this work, but may be of interest for other purposes.

Results in the paper can be re-produced by running
```bash
  sh code/run.sh <CORES> <N_QMC>
```
Where `CORES` is the number of processors you want to use to run the script and `N_QMC` is the number of quasi Monte Carlo samples. We recommend using `CORES=8` and `N_QMC=1000` as an initial simulation. Results in the paper are based on `10 000 000` QMC simulations, we recommend using HPC clusters in this case.

We also provide a pre-build Docker image which can be used to run the the code in this repository. First thing you need to do is to ensure that you have [docker installed](https://docs.docker.com/get-docker/).

### Docker

To start an interactive docker container you can execute the following command

```bash
docker run --rm -it ghcr.io/erasdna/Parameter-estimation-in-the-dura:main
```
This will re-run the segmentation pipeline and parameter estimation, doing `1000` QMC simulations on `8` cores. More cores and QMC samples can be specified.

## Set-by-step pipeline

We provide a bash script `code/run.sh` to run the entire segmentation and parameter estimation pipeline automatically

### Segmentation of the cranial dura

Our segmentation of the cranial dura of Gonzo can be run using the `code/dura.smk` file:

```bash
snakemake -s code/dura.smk --cores 1
```

This will produce segmentations, with corresponding look-up tables (freeview-compatible) and statistics in the `data/paper/dura` folder.

### Parameter estimation

Parameter estimates can be run using:

```bash
snakemake -s code/parameter_estimation.smk --cores 8 --config qmc_samples=1000
```
This will create the reduced model used parameter estimation, run QMC simulations and compute relevant derived quantities.

Relevant outputs can be found in the `data/paper/parameter_estimation` folder.

Please be aware that the results in the full paper are based on 10 million QMC simulations. If you want to run this number of samples, please consider using an HPC cluster resource.


## Citation

```
@software{Girelli2026,
  author = {Girelli, Alberto and Solheim, Andreas and Lysan, Sofie and Storås, Tryggve and Nordengen, Kaja and Mardal, Kent-Andre},
  doi = {10.5281/zenodo.1234},
  month = {5},
  title = {{Multiphasic Transport in the Human Cranial Dura: How Meningeal Lymphatic Vessels Drive Waste Clearance from Cerebrospinal Fluid}},
  url = {https://github.com/Erasdna/Parameter-estimation-in-the-dura},
  version = {1.0},
  year = {2026}
}
```


## Having issues
If you have any troubles please file and issue in the GitHub repository.

## License
MIT
