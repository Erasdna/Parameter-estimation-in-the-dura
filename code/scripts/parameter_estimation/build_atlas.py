import numpy as np
from src.spatial import spherical_patch_mapping, compute_connectivity
from snakemake.script import snakemake
from pathlib import Path
import nibabel as nib
from mritk.segmentation import Segmentation, procedural_freesurfer_lut
from mritk import MRIData
from scipy.ndimage import center_of_mass, sum
import pandas as pd
from scipy.interpolate import NearestNDInterpolator

segmentation_path = Path(snakemake.input.simplified_layer)
segmentation = Segmentation.from_file(segmentation_path)

# Load Data
dura_mask = segmentation.mri.data >= 10
dura_coords = np.argwhere(dura_mask)
csf_mask = (segmentation.mri.data < 10) * (segmentation.mri.data > 0)
csf_coords = np.argwhere(csf_mask)

# Create synthetic spherical atlas
print("Making spherical atlas")
coords_mm = nib.affines.apply_affine(segmentation.mri.affine, dura_coords)
patch_ids = spherical_patch_mapping(coords_mm, n_theta=10, n_phi=15)

synthetic_atlas = np.zeros_like(segmentation.mri.data, dtype=np.int32)
synthetic_atlas[*dura_coords.T] = patch_ids
patches = np.unique(patch_ids)

# Extrapolate patches to CSF
interpolator = NearestNDInterpolator(
    dura_coords,
    synthetic_atlas[*dura_coords.T] + 1000,
)
csf_values = interpolator(*np.argwhere(csf_mask).T)
synthetic_atlas[*csf_coords.T] = csf_values


print("Saving atlas")
unique_patches = np.unique(synthetic_atlas)
lut = procedural_freesurfer_lut(
    unique_patches,
    ["CSF" if x > 1000 else "Dura" for x in unique_patches],
    cmap="viridis",
)
extended_segmentation = Segmentation(
    MRIData(synthetic_atlas, segmentation.mri.affine),
    lut=lut,
)
extended_segmentation.save(Path(snakemake.output.patch_segmentation))


# Use only dura for connectivity
synthetic_atlas_dura = np.where(dura_mask, synthetic_atlas, 0)

print("Computing center of mass")
# Compute center of each patch
centers = np.array(
    center_of_mass(
        np.ones_like(synthetic_atlas_dura),
        labels=synthetic_atlas_dura,
        index=patches,
    ),
)
physical_center = nib.affines.apply_affine(segmentation.mri.affine, centers) / 100

# Compute volume of each patch
dm_scale = np.linalg.det(segmentation.mri.affine[:3, :3] * 1e-2)
volumes = (
    sum(np.ones_like(synthetic_atlas_dura), synthetic_atlas_dura, patches) * dm_scale
)

# Save dataframe
df = pd.DataFrame(physical_center, columns=[f"Centroid_{x}" for x in ["X", "Y", "Z"]])
df["Volume_dm3"] = volumes
df["ROI_ID"] = patches
df = df[
    ["ROI_ID", "Volume_dm3"] + [f"Centroid_{x}" for x in ["X", "Y", "Z"]]
].set_index("ROI_ID")
df.to_csv(Path(snakemake.output.nodes), index=True)

# Make connectivity matrix
print("Computing connectivity matrix")
interface_size = (segmentation.mri.affine[0, 0] / 100) ** 2
edge_df = compute_connectivity(synthetic_atlas_dura, interface_size)
left_centroids = df.loc[edge_df["Source_ROI"].values].values[:, 1:]
right_centroids = df.loc[edge_df["Target_ROI"].values].values[:, 1:]
distance = np.linalg.norm(right_centroids - left_centroids, axis=1)
edge_df["Distance_dm"] = distance

edge_df.to_csv(Path(snakemake.output.edges), index=False)
