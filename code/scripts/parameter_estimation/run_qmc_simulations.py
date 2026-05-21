import os

# Prevent oversubcribing hpc resources
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"

import numpy as np
import pandas as pd
from snakemake.script import snakemake
from src.modeling import (
    fit_CSF_data,
    setup_matrices,
    gamma_func,  # noqa: F401
    double_exp_func,  # noqa: F401
)
import src.constants as constants
from src.qmc import run_qmc
from src.solver import simulate
from functools import partial

NUM_CORES = snakemake.threads

BOUNDS = {
    "phi": snakemake.params.phi,
    "P_bot": snakemake.params.P_bot,
    "D_L": snakemake.params.D_L,
    "P_DcLV": snakemake.params.P_DcLV,
    "phi_L": snakemake.params.phi_L,
}
PARAM_NAMES = list(BOUNDS.keys())
N_QMC = snakemake.params.qmc_samples

# Data reading
df_nodes = pd.read_csv(snakemake.input.nodes)
df_edges = pd.read_csv(snakemake.input.edges)
df_concentration = pd.read_csv(snakemake.input.concentration).rename(
    columns={"label": "ROI_ID"},
)
df_dura = df_concentration.query("ROI_ID<1000")
df_csf = df_concentration.query("ROI_ID>1000")
df_csf["ROI_ID"] -= 1000

N_NODES = len(df_nodes)


def clip_invalid(arr):
    arr = np.where(arr < 0, 0, arr)
    arr = np.where(np.isnan(arr), 0, arr)
    return arr


dura_data = {}
for statistic in ["median", "PC25", "PC75"]:
    data = np.zeros((N_NODES, 5))
    data[:, 1:] = df_dura.query(f"statistic=='{statistic}'").pivot(
        values="value",
        index="ROI_ID",
        columns="session",
    )
    dura_data[statistic] = clip_invalid(data)

csf_data_median = np.zeros((N_NODES, 5))
csf_data_median[:, 1:] = clip_invalid(
    df_csf.query("statistic=='median'")
    .pivot(values="value", index="ROI_ID", columns="session")
    .to_numpy(),
)

# Interpolate CSF in time
fitted_csf_functions = fit_CSF_data(csf_data_median, constants.TIME_POINTS)
# Build interaction matrices
M_G_dense, M_laplacian = setup_matrices(df_edges, df_nodes)

simulator = partial(
    simulate,
    (M_laplacian, M_G_dense),
    dura_data,
    fitted_csf_functions,
    df_nodes["Centroid_Z"].values,
)
df = run_qmc(simulator, NUM_CORES, BOUNDS, N_QMC)
df.to_csv(snakemake.output.simulation_outcomes, index=False)

# --- REPORT IQR ---
total_pts = N_NODES * 4
threshold_90 = total_pts * 0.90

print("\n" + "=" * 50)
print(">>> QMC IQR COVERAGE REPORT")
print("=" * 50)
print(f"Total Valid Runs: {len(df)} / {N_QMC}")
print(f"Max Points in IQR: {df.iloc[0]['Pts_in_IQR']} / {total_pts}")
print(f"  - Runs with 100% points in IQR: {len(df[df['Pts_in_IQR'] == total_pts])}")
print(f"  - Runs with > 90% points in IQR: {len(df[df['Pts_in_IQR'] >= threshold_90])}")
print(
    f"  - Runs with > 80% points in IQR: {len(df[df['Pts_in_IQR'] >= total_pts * 0.8])}",
)
print(
    f"  - Runs with > 50% points in IQR: {len(df[df['Pts_in_IQR'] >= total_pts * 0.5])}",
)
print("=" * 50)
