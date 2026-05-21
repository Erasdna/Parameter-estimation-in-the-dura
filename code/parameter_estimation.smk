configfile: "code/configs/config.yaml"


GONZO = config["data_folder"]  # Path to the gonzo dataset
T1MAPS = GONZO + "T1maps/"

SEGMENTATIONS = config["output_folder"] + "dura/"
OUTPUT_FOLDER = config["output_folder"] + "parameter_estimation/"


rule all:
    input:
        OUTPUT_FOLDER + "patch_segmentation.nii.gz",
        OUTPUT_FOLDER + "concentration.csv",
        OUTPUT_FOLDER + "simulations.csv",
        expand(
            OUTPUT_FOLDER + "derived/{files}.csv",
            files=[
                "top_ensemble_derived_quantities",
                "statistics_derived_quantities",
                "parameter_CIs",
                "parameter_statistics",
            ],
        ),


rule prepare_data:
    input:
        simplified_layer=SEGMENTATIONS + "segmentations/dura_simplified_layer.nii.gz",
    output:
        patch_segmentation=OUTPUT_FOLDER + "patch_segmentation.nii.gz",
        patch_lut=OUTPUT_FOLDER + "patch_segmentation.csv",
        nodes=OUTPUT_FOLDER + "nodes.csv",
        edges=OUTPUT_FOLDER + "edges.csv",
    script:
        "scripts/parameter_estimation/build_atlas.py"


rule compute_concentration:
    input:
        simplified_layer=SEGMENTATIONS + "segmentations/dura_simplified_layer.nii.gz",
        patch_segmentation=OUTPUT_FOLDER + "patch_segmentation.nii.gz",
        patch_lut=OUTPUT_FOLDER + "patch_segmentation.csv",
        T1maps=expand(
            T1MAPS + "sub-01_{session}_T1map_hybrid.nii.gz",
            session=["ses-01", "ses-02", "ses-03", "ses-04", "ses-05"],
        ),
    output:
        stats=OUTPUT_FOLDER + "concentration.csv",
    params:
        T1folder=T1MAPS,
    script:
        "scripts/parameter_estimation/compute_concentration.py"


rule run_simulations:
    input:
        nodes=OUTPUT_FOLDER + "nodes.csv",
        edges=OUTPUT_FOLDER + "edges.csv",
        concentration=OUTPUT_FOLDER + "concentration.csv",
    output:
        simulation_outcomes=OUTPUT_FOLDER + "simulations.csv",
    params:
        phi=config["phi"],
        P_bot=config["P_bot"],
        D_L=config["D_L"],
        P_DcLV=config["P_DcLV"],
        phi_L=config["phi_L"],
        qmc_samples=config["qmc_samples"],
    threads: workflow.cores
    script:
        "scripts/parameter_estimation/run_qmc_simulations.py"


rule compute_dervied_quantities:
    input:
        nodes=OUTPUT_FOLDER + "nodes.csv",
        edges=OUTPUT_FOLDER + "edges.csv",
        simulation_outcomes=OUTPUT_FOLDER + "simulations.csv",
    output:
        expand(
            OUTPUT_FOLDER + "derived/{files}.csv",
            files=[
                "top_ensemble_derived_quantities",
                "statistics_derived_quantities",
                "parameter_CIs",
                "parameter_statistics",
            ],
        ),
    params:
        savepath=OUTPUT_FOLDER + "derived/",
    script:
        "scripts/parameter_estimation/compute_derived_quantities.py"
