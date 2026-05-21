from snakemake.script import snakemake
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import scienceplots  # noqa: F401
import numpy as np
import nibabel as nib
import matplotlib.colors as mcolors

plt.style.use(["science", "no-latex"])
# Load save path
savepath = Path(snakemake.params.savepath)

# Load previously computed stats and LUT
stats = pd.read_csv(snakemake.input.stats)
lut = pd.read_csv(
    Path(snakemake.input.lut),
    sep="\t",
    index_col=0,
    names=["description", "R", "G", "B", "A"],
).rename_axis("label")

color_list = ["C1", "C2", "C3", "C4", "C6", "C7"]
regions = set(lut.index.values)

color_mapping = {
    "cerebellum_brainstem": color_list[0],
    "anterior-inferior": color_list[5],
    "anterior-superior": color_list[3],
    "posterior-inferior": color_list[1],
    "posterior-superior": color_list[2],
}

region_naming = {
    "cerebellum_brainstem": "Cerebellum + brainstem",
    "anterior-inferior": "Anterior-inferior",
    "anterior-superior": "Anterior-superior",
    "posterior-inferior": "Posterior-inferior",
    "posterior-superior": "Posterior-superior",
}

# Plot median and IQR for each region
stat_df = stats.query("data_type=='hybrid'")

x = [4, 24, 48, 72]

# Plot evolution in each region
for region in region_naming.keys():
    region_lut = lut.query(f"description=='{region}' | description=='{region}-dura'")
    rois = region_lut.index.values

    fig, ax = plt.subplots()

    for roi in rois:
        pivot = stat_df.query(f"label=={roi}").pivot(
            index="session",
            columns="statistic",
            values="value",
        )
        tissue = "CSF" if roi < 10 else "Dura"
        ax.plot(
            x,
            pivot["median"]
            / 3.2,  # Convert to mmol/L using relaxivity of 3.2 L/mmol/s,
            label=f"{tissue}",
            marker="o",
            color="k" if tissue == "Dura" else color_mapping[region],
        )
        ax.fill_between(
            x,
            pivot["PC25"] / 3.2,
            pivot["PC75"] / 3.2,
            alpha=0.3,
            color="grey" if tissue == "Dura" else color_mapping[region],
        )
    ax.set_title(region_naming[region])
    ax.set_xticks([4, 24, 48, 72], labels=["4h", "24h", "48h", "72h"])
    leg = ax.legend()
    leg.set_gid("legend")
    ax.set_xlabel("Time after contrast injection")
    ax.set_ylabel("Gd concentration [mmol/L]")
    fig.savefig(savepath.joinpath(f"{region}.png"), dpi=300)
    fig.savefig(savepath.joinpath(f"{region}.pdf"), dpi=300)
    plt.close()

# Plot evolution in CSF vs dura
stat_df = stats.query("statistic=='median' & data_type=='hybrid'")

for tissue_type in ["csf", "dura"]:
    fig, ax = plt.subplots()
    for region in region_naming.keys():
        use_region = region + "-dura" if tissue_type == "dura" else region
        region_lut = lut.query(f"description=='{use_region}'")
        roi = region_lut.index.values[0]

        pivot = stat_df.query(f"label=={roi}").pivot(
            index="session",
            columns="statistic",
            values="value",
        )
        ax.plot(
            x,
            pivot["median"] / 3.2,
            label=f"{region}",
            marker="o",
            color=color_mapping[region],
        )

    ax.set_xticks([4, 24, 48, 72], labels=["4h", "24h", "48h", "72h"])
    ax.set_xlabel("Time after contrast injection")
    ax.set_ylabel("Gd concentration [mmol/L]")
    fig.savefig(savepath.joinpath(f"{tissue_type}.png"), dpi=300)
    fig.savefig(savepath.joinpath(f"{tissue_type}.pdf"), dpi=300)
    ax.legend(
        bbox_to_anchor=(1.05, 1),
        loc="upper left",
        labels=list(region_naming.values()),
    )
    fig.savefig(savepath.joinpath(f"{tissue_type}_label.png"), dpi=300)
    plt.close()

# Plot regions on brain

seg = nib.load(snakemake.input.simple_segmentation).get_fdata().astype(int)
layer_seg = nib.load(snakemake.input.simple_segmentation_layer).get_fdata().astype(int)

fig, ax = plt.subplots(figsize=(10, 10))

seg = np.where(seg >= 10, seg // 10, seg)  # Remove dura/csf distinction for plotting
layer_seg = np.where(
    layer_seg >= 10,
    layer_seg // 10,
    layer_seg,
)

background = np.ones_like(seg) * np.inf

custom_cmap = mcolors.ListedColormap(color_list)
boundaries = [0.5, 1.5, 2.5, 3.5, 4.5, 5.5]
norm = mcolors.BoundaryNorm(boundaries, custom_cmap.N)

ax.imshow(
    np.rot90(background[193]),
    cmap="gray",
    # alpha=0.3,
    gid="background",
    interpolation="nearest",
    rasterized=True,
    vmin=0,
    vmax=1,
)
ax.pcolormesh(
    np.ma.masked_values(np.rot90(seg[193]), 0),
    cmap=custom_cmap,
    norm=norm,
    alpha=0.5,
    gid="segmentation",
    rasterized=True,
)
ax.pcolormesh(
    np.ma.masked_values(np.rot90(layer_seg[193]), 0),
    cmap=custom_cmap,
    norm=norm,
    alpha=1,
    gid="segmentation",
    rasterized=True,
)

ax.axis("off")
fig.savefig(
    fname=savepath.joinpath("regions.png"),
    transparent=False,
    dpi=500,
)
fig.savefig(
    savepath.joinpath("regions.svg"),
    format="svg",
    transparent=True,
)

# Work-around for snakemake
fig.savefig(
    fname=savepath.joinpath("done.png"),
    transparent=False,
    dpi=50,
)
