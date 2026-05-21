import numpy as np
from scipy.ndimage import binary_dilation
import pandas as pd
from tqdm import tqdm


def spherical_patch_mapping(coords_mm, n_theta, n_phi):
    """Maps 3D coordinates to spherical patches."""
    center = np.mean(coords_mm, axis=0)
    rel_coords = coords_mm - center
    r = np.linalg.norm(rel_coords, axis=1)

    # Avoid division by zero
    theta = np.arccos(np.clip(rel_coords[:, 2] / np.where(r == 0, 1e-10, r), -1.0, 1.0))
    phi = np.mod(np.arctan2(rel_coords[:, 1], rel_coords[:, 0]), 2 * np.pi)

    t_idx = np.clip(
        np.digitize(theta, np.linspace(0, np.pi, n_theta + 1)) - 1,
        0,
        n_theta - 1,
    )
    p_idx = np.clip(
        np.digitize(phi, np.linspace(0, 2 * np.pi, n_phi + 1)) - 1,
        0,
        n_phi - 1,
    )

    return t_idx * n_phi + p_idx + 1


def get_bounding_box(mask):
    """Returns the bounding box limits of a 3D boolean mask."""
    coords = np.argwhere(mask > 0)
    coord_min, coord_max = coords.min(axis=0), coords.max(axis=0)
    arr = []
    for i in range(len(coord_min)):
        arr.append(coord_min[i])
        arr.append(coord_max[i])
    return np.array(arr)


def compute_connectivity(label_array, interface_size):
    # Depth interfaces (Z-axis / axis 0)
    # Compare each slice to the slice immediately below it
    z_front = label_array[:-1,]
    z_back = label_array[1:,]
    mask_z = z_front != z_back
    pairs_z = np.column_stack((z_front[mask_z], z_back[mask_z]))

    # Vertical interfaces (Y-axis / axis 1)
    # Compare each row to the row immediately below it
    y_top = label_array[:, :-1]
    y_bottom = label_array[:, 1:]
    mask_y = y_top != y_bottom
    pairs_y = np.column_stack((y_top[mask_y], y_bottom[mask_y]))

    # Horizontal interfaces (X-axis / axis 2)
    # Compare each column to the column immediately to its right
    x_left = label_array[..., :-1]
    x_right = label_array[..., 1:]
    mask_x = x_left != x_right
    pairs_x = np.column_stack((x_left[mask_x], x_right[mask_x]))

    # Combine all touching pairs from all three dimensions
    all_pairs = np.vstack((pairs_z, pairs_y, pairs_x))

    # Sort each pair so that order doesn't matter (A, B) == (B, A)
    all_pairs.sort(axis=1)

    # Filter out the background
    valid_mask = (all_pairs[:, 0] != 0) & (all_pairs[:, 1] != 0)
    patch_pairs = all_pairs[valid_mask]

    # Count the unique pairs to get the interface area (in voxel faces)
    unique_pairs, counts = np.unique(patch_pairs, axis=0, return_counts=True)

    # Remove corners included unintentionally:
    edges_data = {"Source_ROI": [], "Target_ROI": [], "Contact_Area_dm2": []}

    for (roi_A, roi_B), count in tqdm(zip(unique_pairs, counts), total=len(counts)):
        if count <= 1:
            continue
        if len(edges_data["Source_ROI"]) < 1 or roi_A != edges_data["Source_ROI"][-1]:
            roi_A_dilated = binary_dilation(label_array * (label_array == roi_A))

        roi_B_mask = label_array == roi_B
        edges_data["Contact_Area_dm2"].append(
            (roi_A_dilated * roi_B_mask).sum() * interface_size,
        )
        edges_data["Source_ROI"].append(roi_A)
        edges_data["Target_ROI"].append(roi_B)
    return pd.DataFrame(edges_data)
