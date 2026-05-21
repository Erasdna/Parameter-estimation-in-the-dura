<!-- cspell:ignore millimoles mmol numpy nibabel scipy scikit matplotlib Allwmake topo Allrun Allclean intrathecal gadobutrol -->

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
Our analysis suggests that purely diffusive macroscopic models fail to replicate in vivo tracer distribution. Conversely, a dual-porosity advection-diffusion framework successfully captures the tissue kinetics. The transport within the mLV compartment is characterized by highly convective clearance (with $Pe>1$), with macroscopic fluid velocities and effective lymphatic diffusion coefficients confirming an advective capacity sufficient all the total albumin content in the Cerebrospinal Fluid (CSF).

**Conclusion:**
-
Solute clearance within the human cranial dura mater is a multiphasic, advection-dominated process. Rather than relying on passive Fickian diffusion, molecules are actively and rapidly cleared through the meningeal lymphatic network via convective mechanisms. This hierarchical modeling approach bridges clinical imaging and microscopic fluid dynamics, clarifying the primary mechanical drivers of brain homeostasis.


# Reproducing

This repository contains two main parts:
- Under `code/` and `src/` you will find the scripts and code necessary to construct a segmentation of the cranial dura and estimates of various transport properties in the dura. The `docker/` contains a Docker container which can be built to reproduce results
- Under `openFoam_simulation` you will find the code necessary to run full 3D simulations of the cranial dura. Running this code requires the installation of OpenFoam (dev).

## Segmentation and parameter estimation

We provide the tools to reproduce our segmentation of the cranial dura and the associated parameter estimates. Instructions for installing necessary tools are provided, as well as a Docker file to run the entire project.

---

### Docker
The easiest way to reproduce the results in this repo is to use Docker. We provide a pre-built docker container you can use. First thing you need to do is to ensure that you have [docker installed](https://docs.docker.com/get-docker/).

To start an interactive docker container you can execute the following command

```bash
docker run --rm -it ghcr.io/erasdna/duramodels:main
```
This will re-run the segmentation pipeline and parameter estimation, doing `1000` QMC simulations on `8` cores. More cores and QMC samples can be specified.

---

If you do not want to use Docker to re-run the pipeline, we invite you to follow the following setup steps:

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

### Set-by-step pipeline

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



## 3D simulations in openFroam

OpenFoam Case contains the custom OpenFOAM solver (`doubleDarcyStarlingFoam`), the simulation setup, automation scripts, and the Python mapping tool required to run the simulations.

---

### Important Note on Units of Measurement
Please pay close attention when interpreting the simulation results, editing dictionaries, or reviewing the code. The units used in this OpenFOAM case are **not** the standard SI units. Instead, the following unit system is applied:
* **Length:** millimeters (mm)
* **Time:** hours (h)
* **Mass:** grams (g)
* **Amount of substance:** millimoles (mmol)

*Why?* This custom unit system was adopted strictly for convenience to perfectly match the scale of the provided clinical MRI data and the imported STL geometries, avoiding rounding errors during conversions.

---

## Pre-meshed and Compiled Dataset (Zenodo)

If you wish to skip the meshing process and the boundary condition mapping, the fully meshed case,along with the pre-mapped boundary conditions and the solutions, are freely available on Zenodo.

**Download the complete dataset here:** https://doi.org/10.5281/zenodo.20120581

Or run the following bash command to automatically download and move all the mesh and BC in the right folder for the case "Data-integrated 72 hour simulation of intrathecal gadobutrol"

```bash
  sh openFoam_simulation/run_Mapping3D.sh
```

Or run the following bash command to automatically download and move all the mesh and BC in the right folder for the case "Clearance of albumin through dural lymphatics"

```bash
  sh openFoam_simulation/run_Albumin3D.sh
```

If you want to Run the meshing pipes, please run the following command to directly download and move .stl files and initial/boundary conditions folder, specifying which case you are interested in (Albumin/Mapping):
```bash
  sh openFoam_simulation/run.sh <Albumin/Mapping>
```

---

## How to Run the Case from Scratch

If you prefer to generate the setup from the source files provided in this repository, please follow the steps below. If not, skip directly to the running case. Please ensure that [OpenFoam](https://openfoam.org/version/dev/) is installed on your system and sourced.

### Compile the Custom Solver
Before running the simulation, you must compile the custom solver `doubleDarcyStarlingFoam`:
```bash
./Allwmake
```

### Mesh Generation
Navigate to the root folder of the simulation case and execute the automated meshing script. This will run `surfaceFeatureExtract`, `pMesh`, `topoSet`, `setFields`, and `createPatch`.
```bash
cd Dura_simulations
./Allrun.mesh
```

### Mapping MRI Data to OpenFOAM
To generate the time-varying boundary conditions (Cin_field) from your clinical NIfTI images, you must use the provided Python script.

⚠️ Attention: Before running the script, open map_MRI_to_openFoam.py with a text editor and ensure that the base directory variable points exactly to the folder where your .nii.gz NIfTI files and segmentations are saved on your local machine.

Once configured, run the script:
```bash
uv run python map_MRI_to_openFoam.py
```

### Run Simulations
Once the mesh is generated and the boundary conditions are mapped, you can start the parallel simulation. The script will automatically decompose the domain, run the solver in parallel, and reconstruct the results:
```bash
./Allrun.solve
```
⚠️ Attention: We suggest to run the simulation in an HPC cluster resource.

### Clean the Directory
To reset the case to its initial state and remove all generated mesh, processor folders, and log files, run:
```bash
./Allclean
```


## Citation

```
@software{Girelli2026,
  author = {Girelli, Alberto and Solheim, Andreas and Lysan, Sofie and Storås, Tryggve and Nordengen, Kaja and Mardal, Kent-Andre},
  doi = {10.5281/zenodo.1234},
  month = {5},
  title = {{Multiphasic Transport in the Human Cranial Dura: How Meningeal Lymphatic Vessels Drive Waste Clearance from Cerebrospinal Fluid}},
  url = {https://github.com/scientificcomputing/example-paper},
  version = {1.0},
  year = {2026}
}
```


## Having issues
If you have any troubles please file and issue in the GitHub repository.

## License
MIT
