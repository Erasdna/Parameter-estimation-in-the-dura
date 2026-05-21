from snakemake.script import snakemake
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import scienceplots  # noqa: F401
import numpy as np
from mritk.r1 import convert_t1_to_r1
from mritk import MRIData
from mritk.segmentation import Segmentation
import seaborn as sns
from scipy.stats import wilcoxon, quantile_test

plt.style.use(["science", "no-latex"])

filename = lambda ses, data_type: f"sub-01_{ses}_T1map_{data_type}.nii.gz"

# Load previously computed stats and LUT
lut = pd.read_csv(
    Path(snakemake.input.lut),
    sep="\t",
    index_col=0,
    names=["description", "R", "G", "B", "A"],
).rename_axis("label")
seg_data = MRIData.from_file(Path(snakemake.input.segmentation), orient=False)
seg = Segmentation(seg_data, lut=lut)

color_list = ["C1", "C2", "C3", "C4", "C6", "C7"]
regions = set(lut.index.values)

T1folder = Path(snakemake.params.T1folder)
T1types = snakemake.params.T1map_types

sequence_naming = {"looklocker": "Look-Locker", "mixed": "Mixed", "hybrid": "Hybrid"}

final_df: pd.DataFrame = None
for session in ["ses-02", "ses-03", "ses-04", "ses-05"]:
    for data_type in T1types:
        # Find ses-01 data
        ses01_path = T1folder.joinpath(filename("ses-01", data_type))
        ses01 = MRIData.from_file(ses01_path, orient=False)
        baseline_R1 = convert_t1_to_r1(ses01)

        print(f"Collecting data for session {session}, {data_type}")
        ses = MRIData.from_file(
            T1folder.joinpath(filename(session, data_type)),
            orient=False,
        )
        Gd_concentration = convert_t1_to_r1(ses).data - baseline_R1.data
        Gd_concentration = np.where(np.isnan(Gd_concentration), 0, Gd_concentration)
        csf_data = Gd_concentration[(seg.mri.data < 10) * (seg.mri.data > 0)]
        csf_df = pd.DataFrame(csf_data, columns=["value"])
        csf_df["tissue"] = "CSF"

        dura_data = Gd_concentration[seg.mri.data >= 10]
        dura_df = pd.DataFrame(dura_data, columns=["value"])
        dura_df["tissue"] = "Dura"

        df = pd.concat([dura_df, csf_df], ignore_index=True)
        df["data_type"] = sequence_naming[data_type]
        df["session"] = session

        if final_df is None:
            final_df = df.copy()
        else:
            final_df = pd.concat([final_df, df], ignore_index=True)
        print(final_df)

        for data, name in zip([dura_data, csf_data], ["dura", "csf"]):
            wilcoxon_stat, wilcoxon_p = wilcoxon(
                data,
                alternative="greater",
            )
            q_test = quantile_test(data, alternative="greater")
            print(f"Wilcoxon signed-rank test in the {name}: p={wilcoxon_p}")
            print(f"Quantile test in the {name}: p={q_test.pvalue}")

cm = 1 / 2.54

final_df["value"] = final_df["value"] / 3.2
for i, tissue in enumerate(["CSF", "Dura"]):
    fig, ax = plt.subplots(figsize=(17 * cm, 7 * cm))

    boxplot = sns.boxplot(
        data=final_df.query(f"tissue=='{tissue}'"),
        x="session",
        y="value",
        hue="data_type",
        palette="tab10",
        ax=ax,
        showfliers=False,
        gap=0.1,
        legend="auto" if i == 0 else False,
    )

    ax.set_ylabel("Gd concentration [mmol/L]")
    if i == 0:
        ax.legend(title=r"$T_1$ sequence")
    ax.set_xticks([0, 1, 2, 3], [r"4h", r"24h", r"48h", r"72h"])

    ax.set_xlabel("Time after injection")

    fig.savefig(
        Path(snakemake.params.savepath).joinpath(tissue).with_suffix(".png"),
        dpi=300,
    )

    fig.savefig(
        Path(snakemake.params.savepath).joinpath("done").with_suffix(".png"),
        dpi=50,
    )

    plt.close()
