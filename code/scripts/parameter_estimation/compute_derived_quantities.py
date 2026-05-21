import pandas as pd
from snakemake.script import snakemake
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import scipy.sparse.linalg as spla
import src.constants as constants
import scipy.sparse as sp
from src.modeling import setup_matrices

# Read simulation outcomes
df = pd.read_csv(Path(snakemake.input.simulation_outcomes), index_col=0)
PARAM_NAMES = df.columns[:5]

# Save path
OUTPUT_DIR = Path(snakemake.params.savepath)

df_nodes = pd.read_csv(snakemake.input.nodes)
df_edges = pd.read_csv(snakemake.input.edges)

# We run statistical analysis of the 200 top parameter fits
top_ensemble = df.sort_values("Weighted_Loss", ascending=True).head(200)

# new stats
print("\n>>> Stats Description on the Top Fits...")
# all principal stats
stats_df = (
    top_ensemble[PARAM_NAMES].describe(percentiles=[0.025, 0.25, 0.5, 0.75, 0.975]).T
)

stats_df.rename(
    columns={
        "mean": "Mean",
        "std": "Std_Dev",
        "min": "Min",
        "50%": "Median",
        "max": "Max",
    },
    inplace=True,
)

# save basic statistics on parameters in top ensemble
stats_df["IQR"] = stats_df["75%"] - stats_df["25%"]
stats_df.to_csv(
    OUTPUT_DIR.joinpath("parameter_statistics.csv"),
    index_label="Parameter",
)

# Save confidence intervals for parameters in top ensemble
ci = {
    n: [top_ensemble[n].quantile(0.025), top_ensemble[n].quantile(0.975)]
    for n in PARAM_NAMES
}
pd.DataFrame(ci, index=["Lower_95", "Upper_95"]).to_csv(
    OUTPUT_DIR.joinpath("parameter_CIs.csv"),
)

# compute derived quantities from estimated parameters
Z_coords = df_nodes["Centroid_Z"].values
N_NODES = len(Z_coords)
node_top_idx, node_bot_idx = np.argmax(Z_coords), np.argmin(Z_coords)
_, M_laplacian = setup_matrices(df_edges, df_nodes)

L_head = Z_coords.max() - Z_coords.min()


# Base pressure field
A = constants.K_DARCY_FIXED * M_laplacian - sp.eye(N_NODES) * constants.GAMMA_FIXED
b = -constants.GAMMA_FIXED * constants.H_D_FIXED * np.ones(N_NODES)
A_lil = A.tolil()
A_lil[node_bot_idx, :] = 0
A_lil[node_bot_idx, node_bot_idx] = 1
b[node_bot_idx] = 0
HL = spla.spsolve(A_lil.tocsr(), b)
Actual_H_max = np.max(HL)
L_head = Z_coords.max() - Z_coords.min()

current_u_interst = (constants.K_DARCY_FIXED * Actual_H_max / L_head) / top_ensemble[
    "phi_L"
]

derived_metrics = {
    "Idx": top_ensemble.index.values,
    "Lymphatic_Peclet": (current_u_interst * L_head) / top_ensemble["D_L"],
    "Biot_Dura": (top_ensemble["P_bot"] * 0.01) / constants.D_D,
    "Advection_Time_h": np.where(
        current_u_interst > 0,
        L_head / current_u_interst,
        np.inf,
    ),
    "Diffusion_Time_h": (L_head**2) / top_ensemble["D_L"],
    "Weighted_Loss": top_ensemble["Weighted_Loss"],
}

# DataFrame with all the rows
df_derived = pd.DataFrame(derived_metrics)

# Stats (Mean, Std, Min, Max, Quartiles)
# Exclude Run_Idx e Weighted_Loss
df_derived_stats = df_derived.drop(columns=["Idx", "Weighted_Loss"]).describe().T

# Add IQR column
df_derived_stats["IQR"] = df_derived_stats["75%"] - df_derived_stats["25%"]
df_derived_stats.insert(0, "Physical_Quantity", df_derived_stats.index)

# Save
df_derived_stats.to_csv(
    OUTPUT_DIR.joinpath("statistics_derived_quantities.csv"),
    index=False,
)
df_derived.to_csv(
    OUTPUT_DIR.joinpath("top_ensemble_derived_quantities.csv"),
    index=False,
)

print("\n>>> Derived Physical Stats (Péclet, Biot, Tempi) computed and saved on CSV!")

print("\n>>> Generazione Log-Correlation Matrix e PairPlot...")
# copy to not overwrite original data
log_ensemble = top_ensemble[PARAM_NAMES].copy()

for col in PARAM_NAMES:
    if log_ensemble[col].max() / np.maximum(log_ensemble[col].min(), 1e-9) > 100.0:
        log_ensemble[col] = np.log10(log_ensemble[col])
        log_ensemble.rename(columns={col: f"log10_{col}"}, inplace=True)

# log corr matrix (better in this case)
fig, ax = plt.subplots(figsize=(12, 10))
sns.heatmap(
    log_ensemble.corr(),
    annot=True,
    cmap="coolwarm",
    fmt=".2f",
    vmin=-1,
    vmax=1,
)
plt.title(f"Log-Parameter Correlation Matrix (Top {len(top_ensemble)} Runs)")
plt.tight_layout()
plt.savefig(OUTPUT_DIR.joinpath("Correlation_Matrix_Log_QMC.png"), dpi=300)
plt.close()

# PAIRPLOT
sns.pairplot(
    log_ensemble,
    kind="scatter",
    diag_kind="kde",
    plot_kws={"alpha": 0.4, "s": 15, "edgecolor": "none", "color": "red"},
    diag_kws={"fill": True, "color": "red"},
)
fig.suptitle("Parameter Space Exploration (PairPlot)", y=1.02, fontsize=16)
plt.savefig(OUTPUT_DIR.joinpath("Parameter_PairPlot_QMC.png"), dpi=300)
plt.close()
