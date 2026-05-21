# ruff: noqa: E501
import pandas as pd
import numpy as np
import os
import nibabel as nib
import scipy.ndimage as ndi
from sklearn.neighbors import NearestNeighbors
import matplotlib.pyplot as plt
import warnings

warnings.filterwarnings("ignore", category=RuntimeWarning)

# ==============================================================================
# CONFIGURATION
# ==============================================================================
# Define patches and time steps to process
patch_names = ["internalPatch", "AG"]
time_steps = ["0", "4", "24", "48", "72"]

# Map time steps to their respective NIfTI files
time_to_nifti = {
    "4": "ses-02_looklocker_R1.nii.gz",
    "24": "ses-03_looklocker_R1.nii.gz",
    "48": "ses-04_looklocker_R1.nii.gz",
    "72": "ses-05_looklocker_R1.nii.gz",
}

# Base directory for MRI data
base_dir_nifti = "./MRI_data"

# Principal parameters
FS_TAG_MIN = 10000
FS_TAG_MAX = 20000
T1_THRESHOLD_MS = 1500
MAX_DIST_CSF_TO_DURA_MM = 2.0
MAPPING_RADIUS_MM = 0.2  # Search radius for CSF voxels around mesh points
N_SMOOTH_NEIGHBORS = 0  # Number of nearest neighbors for smoothing (0 = disabled)

# ==============================================================================
# READ OPENFOAM MESH AND APPLY INVERSE TRANSFORMATION (OpenFOAM -> MRI Space)
# ==============================================================================
print("Reading OpenFOAM face centers from cell_centers.csv...")
df_centers = pd.read_csv("cell_centers.csv")

openfoam_points = df_centers[["Points:0", "Points:1", "Points:2"]].values
num_points = len(openfoam_points)

print("Applying inverse transformation (OpenFOAM -> MRI space)...")
angle_degrees = -94.3369
offset_array = np.array([3.3959, -0.5899, 19.6657])

theta = np.radians(angle_degrees)
c, s = np.cos(theta), np.sin(theta)
R = np.array(
    [
        [1, 0, 0],
        [0, c, -s],
        [0, s, c],
    ],
)

mri_points = (openfoam_points - offset_array).dot(R)

# ==============================================================================
# CREATE CSF MASK NEAR THE DURA MATER
# ==============================================================================
print(f"\nCreating CSF mask within {MAX_DIST_CSF_TO_DURA_MM} mm of the dura...")
mask_path_dura = os.path.join(
    base_dir_nifti,
    "segmentations/dura_files/dura_mask.nii.gz",
)
mask_path_csf = os.path.join(base_dir_nifti, "segmentations/dura_files/CSF_mask.nii.gz")
ref_nifti_path = os.path.join(base_dir_nifti, "concentration", time_to_nifti["4"])

ref_img = nib.load(ref_nifti_path)
reference_affine = ref_img.affine
voxel_sizes = ref_img.header.get_zooms()[:3]

# Compute distance from the dura mask
mask_dura_data = nib.load(mask_path_dura).get_fdata()
plot_dura_mask = mask_dura_data > 0
dist_map_to_dura = ndi.distance_transform_edt(~plot_dura_mask, sampling=voxel_sizes)

# Filter FreeSurfer CSF mask by tags
mask_csf_data = nib.load(mask_path_csf).get_fdata()
plot_csf_mask_base = (mask_csf_data >= FS_TAG_MIN) & (mask_csf_data <= FS_TAG_MAX)

# Final CSF mask: correct tags, within distance threshold, excluding dura voxels
final_csf_mask = (
    plot_csf_mask_base
    & (dist_map_to_dura <= MAX_DIST_CSF_TO_DURA_MM)
    & (~plot_dura_mask)
)

coords_voxel_csf = np.column_stack(np.where(final_csf_mask))
coords_physical_csf = nib.affines.apply_affine(reference_affine, coords_voxel_csf)

print(f" -> Valid CSF voxels found: {len(coords_voxel_csf)}")

if len(coords_voxel_csf) == 0:
    raise RuntimeError(
        "No CSF voxels found! Check your masks or increase MAX_DIST_CSF_TO_DURA_MM.",
    )

# ==============================================================================
# KDTREE MAPPING (MRI VOXELS -> MESH POINTS)
# ==============================================================================
print(f"\nBuilding KDTree for CSF voxels (Search radius = {MAPPING_RADIUS_MM} mm)...")
nn_mapper = NearestNeighbors(n_jobs=-1)
nn_mapper.fit(coords_physical_csf)

# Find CSF voxels within the specified radius for each mesh point
print("Mapping voxels to mesh points...")
indices_radius = nn_mapper.radius_neighbors(
    mri_points,
    radius=MAPPING_RADIUS_MM,
    return_distance=False,
)

# Fallback: nearest neighbor for points with no voxels in radius
nn_fallback = NearestNeighbors(n_neighbors=1, n_jobs=-1)
nn_fallback.fit(coords_physical_csf)
_, indices_fallback = nn_fallback.kneighbors(mri_points)

n_fallback = sum(1 for idx in indices_radius if len(idx) == 0)
print(f" -> Points with >=1 voxel in radius: {num_points - n_fallback}/{num_points}")
print(f" -> Points using nearest-neighbor fallback: {n_fallback}")


def median_with_fallback(values_in_radius, fallback_idx, all_values):
    """Returns median of values in radius, or nearest neighbor if radius is empty."""
    if len(values_in_radius) > 0:
        return np.median(values_in_radius)
    else:
        return all_values[fallback_idx[0]]


# ==============================================================================
# RECONSTRUCT CONCENTRATION FIELDS
# ==============================================================================
reconstructed_concs = {}

for t_str, nifti_filename in time_to_nifti.items():
    print(f"\n--- Processing time step: {t_str} h ---")
    ll_path = os.path.join(base_dir_nifti, "concentration", nifti_filename)
    mix_path = os.path.join(
        base_dir_nifti,
        "concentration",
        nifti_filename.replace("looklocker", "mixed"),
    )

    if not os.path.exists(ll_path) or not os.path.exists(mix_path):
        raise FileNotFoundError(f"Missing NIfTI files for time {t_str}.")

    # Select LookLocker or Mixed R1 values based on T1 threshold
    r1_ll = nib.load(ll_path).get_fdata()[final_csf_mask]
    r1_mix = nib.load(mix_path).get_fdata()[final_csf_mask]

    t1_ll = np.divide(1000, r1_ll, out=np.zeros_like(r1_ll), where=r1_ll != 0)
    selected_r1 = np.where(t1_ll < T1_THRESHOLD_MS, r1_ll, r1_mix)
    conc_voxels = selected_r1 * (1.0 / 3.2)  # R1 to concentration conversion

    # Calculate local median for each mesh point
    print(" -> Calculating local median for mesh points...")
    mapped = np.array(
        [
            median_with_fallback(conc_voxels[idx_r], indices_fallback[i], conc_voxels)
            for i, idx_r in enumerate(indices_radius)
        ],
    )

    reconstructed_concs[t_str] = mapped
    print(
        f" -> Stats: min={mapped.min():.4f}, max={mapped.max():.4f}, mean={mapped.mean():.4f}",
    )

# ==============================================================================
# OPTIONAL POINT-TO-POINT SMOOTHING ON OPENFOAM MESH
# ==============================================================================
smoothed_fields = {}

if N_SMOOTH_NEIGHBORS > 1:
    print(f"\nBuilding KDTree for mesh smoothing (K={N_SMOOTH_NEIGHBORS})...")
    nn_smooth = NearestNeighbors(n_neighbors=N_SMOOTH_NEIGHBORS, n_jobs=-1)
    nn_smooth.fit(openfoam_points)
    _, indices_smooth = nn_smooth.kneighbors(openfoam_points)
else:
    indices_smooth = None
    print("\nMesh smoothing disabled.")

for t_str in time_steps:
    print(f"Finalizing concentration field for t={t_str}...")

    if t_str == "0":
        smoothed_fields[t_str] = np.zeros(num_points)
    else:
        raw_conc = reconstructed_concs[t_str] * 1.0e-6  # Unit conversion
        if indices_smooth is not None:
            smoothed_fields[t_str] = np.median(raw_conc[indices_smooth], axis=1)
        else:
            smoothed_fields[t_str] = raw_conc.copy()

    # Export debug points
    pd.DataFrame(
        {
            "x": openfoam_points[:, 0],
            "y": openfoam_points[:, 1],
            "z": openfoam_points[:, 2],
            "concentration": smoothed_fields[t_str],
        },
    ).to_csv(f"DEBUG_points_t{t_str}.csv", index=False)

# ==============================================================================
# WRITE OPENFOAM BOUNDARY DATA
# ==============================================================================
FOAM_HEADER = """\
/*--------------------------------*- C++ -*----------------------------------*\\
  =========                 |
  \\\\      /  F ield         | OpenFOAM
   \\\\    /   O peration     |
    \\\\  /    A nd           |
     \\\\/     M anipulation  |
\\*---------------------------------------------------------------------------*/
"""

for patch in patch_names:
    print(f"\n--- Writing boundaryData for patch: {patch} ---")
    base_dir = f"constant/boundaryData/{patch}"
    os.makedirs(base_dir, exist_ok=True)

    for t_str in time_steps:
        os.makedirs(os.path.join(base_dir, t_str), exist_ok=True)

    # Write points file
    with open(os.path.join(base_dir, "points"), "w") as f:
        f.write(FOAM_HEADER)
        f.write(
            f"FoamFile\n{{\n    version     2.0;\n    format      ascii;\n"
            f'    class       vectorField;\n    location    "constant/boundaryData/{patch}";\n'
            f"    object      points;\n}}\n// * * * //\n\n{num_points}\n(\n",
        )
        for pt in openfoam_points:
            f.write(f"({pt[0]:.6f} {pt[1]:.6f} {pt[2]:.6f})\n")
        f.write(")\n")

    # Write Cin_field file for each time step
    for t_str, conc_values in smoothed_fields.items():
        with open(os.path.join(base_dir, t_str, "Cin_field"), "w") as f:
            f.write(FOAM_HEADER)
            f.write(
                f"FoamFile\n{{\n    version     2.0;\n    format      ascii;\n"
                f'    class       scalarField;\n    location    "constant/boundaryData/{patch}/{t_str}";\n'
                f"    object      Cin_field;\n}}\n// * * * //\n\n{num_points}\n(\n",
            )
            for val in conc_values:
                f.write(f"{val:.6e}\n")
            f.write(")\n")

print("\nSuccess! OpenFOAM boundaryData generated.")

# ==============================================================================
# 3D PREVIEW PLOT
# ==============================================================================
print("\nGenerating 3D preview plot (t=4h)...")
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection="3d")

conc_plot = smoothed_fields["4"]
skip = max(1, len(openfoam_points) // 5000)

sc = ax.scatter(
    openfoam_points[::skip, 0],
    openfoam_points[::skip, 1],
    openfoam_points[::skip, 2],
    c=conc_plot[::skip],
    cmap="viridis",
    s=5,
    alpha=0.8,
)
plt.colorbar(sc, label="CSF Concentration")
ax.set_title("BC Preview t=4h")
ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("Z")

print("Close the plot window to exit.")
plt.show()
