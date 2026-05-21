from snakemake.script import snakemake
from mritk.r1 import convert_t1_to_r1
from mritk import MRIData
from mritk.segmentation import Segmentation
from mritk.statistics.compute_stats import generate_stats_dataframe_rois
from mritk.statistics.stat_functions import PC25, PC75, Median
from pathlib import Path
import pandas as pd
from scipy.stats import wilcoxon, ttest_rel, quantile_test
import numpy as np

T1folder = Path(snakemake.params.T1folder)
T1types = snakemake.params.T1map_types

filename = lambda ses, data_type: f"sub-01_{ses}_T1map_{data_type}.nii.gz"

segmentation_lut = pd.read_csv(
    Path(snakemake.input.lut),
    sep="\t",
    index_col=0,
    names=["description", "R", "G", "B", "A"],
).rename_axis("label")

seg_data = MRIData.from_file(Path(snakemake.input.segmentation), orient=False)
seg = Segmentation(seg_data, lut=segmentation_lut)

final_df = None
for data_type in T1types:
    # Find ses-01 data
    ses01_path = T1folder.joinpath(filename("ses-01", data_type))
    ses01 = MRIData.from_file(ses01_path, orient=False)
    baseline_R1 = convert_t1_to_r1(ses01)

    for session in ["ses-02", "ses-03", "ses-04", "ses-05"]:
        print(f"Computing stats for session {session}, {data_type}")
        ses = MRIData.from_file(
            T1folder.joinpath(filename(session, data_type)),
            orient=False,
        )
        ses_r1 = convert_t1_to_r1(ses)
        Gd_concentration = ses_r1.data - baseline_R1.data
        df = generate_stats_dataframe_rois(
            seg,
            MRIData(Gd_concentration, baseline_R1.affine),
            qois=[PC25, PC75, Median],
        )
        df["session"] = session
        df["data_type"] = data_type

        if final_df is not None:
            final_df = pd.concat([final_df, df], ignore_index=True)
        else:
            final_df = df.copy()

        final_df.to_csv(Path(snakemake.output.stats))

        # Check that signal in the dura is statistically significant:
        print(f"Comparing dura signal to baseline for session {session}")
        rois = np.unique(seg.mri.data[seg.mri.data > 0])
        for roi in rois[rois >= 10]:
            ref = baseline_R1.data[*np.argwhere(seg.mri.data == roi).T]
            session_data = ses_r1.data[*np.argwhere(seg.mri.data == roi).T]
            wilcoxon_stat, wilcoxon_p = wilcoxon(
                session_data,
                ref,
                alternative="greater",
            )
            ttest_rel_stat, ttest_rel_p = ttest_rel(
                session_data,
                ref,
                alternative="two-sided",
            )
            q_test = quantile_test(session_data - ref, alternative="greater")
            print(
                f"Wilcoxon signed-rank test for {roi} in session {session}: W={wilcoxon_stat}, p={wilcoxon_p}, median={np.median(session_data)}",  # noqa: E501
            )
            print(
                f"Paired t-test for {roi} in session {session}: t={ttest_rel_stat}, p={ttest_rel_p}, mean difference={np.mean(session_data)}",  # noqa: E501
            )
            print(
                f"Quantile-test for {roi} in session {session}: p={q_test.pvalue}",
            )
