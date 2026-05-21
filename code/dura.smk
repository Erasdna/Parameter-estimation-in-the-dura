configfile: "code/configs/config.yaml"


GONZO = config["data_folder"]  # Path to the gonzo dataset
OUTPUT_FOLDER = config["output_folder"] + "dura/"
SEGMENTATIONS = GONZO + "segmentations/"
REGISTERED = GONZO + "registered/"
T1MAPS = GONZO + "T1maps/"

data_types = ["looklocker", "mixed"]
T1map_types = data_types + ["hybrid"]

import sys
print("EXECUTABLE:", sys.executable)

rule all:
    input:
        # Combined segmentation
        expand(
            OUTPUT_FOLDER + "segmentations/{files}{suffix}",
            files=[
                "wmparc_tissue_csf_dura",
                "dura_simplified",
                "dura_simplified_layer",
            ],
            suffix=[".nii.gz", ".csv"],
        ),
        OUTPUT_FOLDER + "statistics/stats.csv",
        OUTPUT_FOLDER + "statistics/plots/done.png",
        OUTPUT_FOLDER + "statistics/plots/distributions/done.png",
        expand(
            T1MAPS + "sub-01_{ses}_T1map_{data_type}.nii.gz",
            ses=["ses-01", "ses-02", "ses-03", "ses-04", "ses-05"],
            data_type=["looklocker", "mixed"],
        )

rule make_dura_segmentation:
    input:
        wmparc=SEGMENTATIONS + "sub-01_seg-wmparc_refined.nii.gz",  # wmparc segmentation,
        csf=SEGMENTATIONS + "sub-01_seg-csf-wmparc.nii.gz",  # csf segmentation,
        mixed=REGISTERED + "sub-01_ses-01_acq-mixed_T1map_registered.nii.gz",  # mixed T1map ses 01,
        looklocker=REGISTERED + "sub-01_ses-01_acq-looklocker_T1map_registered.nii.gz",  # looklocker T1map ses 01,
    output:
        expand(
            OUTPUT_FOLDER + "segmentations/{files}{suffix}",
            files=[
                "wmparc_tissue_csf_dura",
                "dura_simplified",
                "dura_simplified_layer",
            ],
            suffix=[".nii.gz", ".csv"],
        ),
    params:
        dura_segmentation=OUTPUT_FOLDER + "segmentations/wmparc_tissue_csf_dura.nii.gz",  # segmentation with parenchyma, csf, and dura,
        simplified_segmentation=OUTPUT_FOLDER + "segmentations/dura_simplified.nii.gz",  # segmentation split into 5 regions, used for plotting,
        simplified_segmentation_layer=OUTPUT_FOLDER
        + "segmentations/dura_simplified_layer.nii.gz",
        # Segmentation including only outer dura and adjacent CSF
    script:
        "scripts/segmentation/make_dura_segmentation.py"


rule copy_T1maps:
    input:
        IN=REGISTERED + "sub-01_{ses}_acq-{data_type}_T1map_registered.nii.gz"
    output:
        OUT=T1MAPS + "sub-01_{ses}_T1map_{data_type}.nii.gz"
    shell:
        "rsync -av {input.IN} {output.OUT}"

rule compute_statistics:
    input:
        T1maps=expand(
            T1MAPS + "sub-01_{session}_T1map_{T1maps}.nii.gz",
            session=["ses-01", "ses-02", "ses-03", "ses-04", "ses-05"],
            T1maps=["looklocker", "mixed", "hybrid"],
        ),
        segmentation=OUTPUT_FOLDER + "segmentations/dura_simplified_layer.nii.gz",
        lut=OUTPUT_FOLDER + "segmentations/dura_simplified_layer.csv",
    output:
        stats=OUTPUT_FOLDER + "statistics/stats.csv",  # Median + IQR stats for dura + csf on the edge
    params:
        T1folder=T1MAPS,
        T1map_types=["looklocker", "mixed", "hybrid"],
    script:
        "scripts/segmentation/compute_statistics.py"


rule plot_statistics:
    input:
        simple_segmentation=OUTPUT_FOLDER + "segmentations/dura_simplified.nii.gz",
        simple_segmentation_layer=OUTPUT_FOLDER
        + "segmentations/dura_simplified_layer.nii.gz",
        stats=OUTPUT_FOLDER + "statistics/stats.csv",
        lut=OUTPUT_FOLDER + "segmentations/dura_simplified_layer.csv",
    output:
        OUTPUT_FOLDER + "statistics/plots/done.png",
    params:
        savepath=OUTPUT_FOLDER + "statistics/plots/",
    script:
        "scripts/segmentation/plot_statistics.py"


rule make_distributions:
    input:
        segmentation=OUTPUT_FOLDER + "segmentations/dura_simplified_layer.nii.gz",
        T1maps=expand(
            T1MAPS + "sub-01_{session}_T1map_{T1maps}.nii.gz",
            session=["ses-01", "ses-02", "ses-03", "ses-04", "ses-05"],
            T1maps=["looklocker", "mixed", "hybrid"],
        ),
        lut=OUTPUT_FOLDER + "segmentations/dura_simplified_layer.csv",
    output:
        OUTPUT_FOLDER + "statistics/plots/distributions/done.png",
    params:
        T1folder=T1MAPS,
        T1map_types=["looklocker", "mixed", "hybrid"],
        savepath=OUTPUT_FOLDER + "statistics/plots/distributions",
    script:
        "scripts/segmentation/plot_distributions.py"
