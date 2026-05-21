from snakemake.script import snakemake
from mritk.segmentation import Segmentation
from mritk import MRIData
from mritk.statistics.compute_stats import generate_stats_dataframe_rois
from mritk.statistics.stat_functions import PC25, PC75, Median
import pandas as pd
from pathlib import Path
from mritk.r1 import convert_t1_to_r1

segmentation_path = Path(snakemake.input.simplified_layer)
segmentation = Segmentation.from_file(segmentation_path)


# Load patches
patch_path = Path(snakemake.input.patch_segmentation)
patch_lut = pd.read_csv(
    Path(snakemake.input.patch_lut),
    sep="\t",
    index_col=0,
    names=["description", "R", "G", "B", "A"],
).rename_axis("label")
patch_data = MRIData.from_file(patch_path, orient=False)
patch_segmentation = Segmentation(patch_data, lut=patch_lut)

filename = lambda ses, data_type: f"sub-01_{ses}_T1map_{data_type}.nii.gz"


# Find ses-01 data
T1folder = Path(snakemake.params.T1folder)
ses01_path = T1folder.joinpath(filename("ses-01", "hybrid"))
ses01 = MRIData.from_file(ses01_path, orient=False)
baseline_R1 = convert_t1_to_r1(ses01)

save_df = None
for session in ["ses-02", "ses-03", "ses-04", "ses-05"]:
    print(f"Computing stats for session {session}")
    ses = MRIData.from_file(
        T1folder.joinpath(filename(session, "hybrid")),
        orient=False,
    )
    ses_r1 = convert_t1_to_r1(ses)
    Gd_concentration = ses_r1.data - baseline_R1.data

    df = generate_stats_dataframe_rois(
        patch_segmentation,
        MRIData(Gd_concentration, baseline_R1.affine),
        qois=[PC25, PC75, Median],
    )
    df["value"] = df["value"] / 3.2
    df["session"] = session

    if save_df is None:
        save_df = df
    else:
        save_df = pd.concat([save_df, df])

    save_df.to_csv(Path(snakemake.output.stats))
