from mritk.segmentation import (
    ExtendedFreeSurferSegmentation,
)
from mritk import MRIData
import numpy as np
from scipy.ndimage import binary_dilation, binary_erosion
from skimage import measure
from scipy.interpolate import NearestNDInterpolator


# Add together to create combined segmentation
def combine_segmentations(seg1, seg2, offset=1000, priority: str = "right"):
    """Combine two segmentations. Labels from seg2 are offset by offset.
    In areas with overlapping labels priority decides which array is used

    Args:
        seg1 (FreeSurferSegmentation): First segmentation
        seg2 (FreeSurferSegmentation): Second segmentation
        offset (int, optional): Offset to be applied to seg2 labels. Defaults to 1000.
        priority (str, optional): Label file which will be given priority on overlaps.
        Defaults "right".

    Returns:
        ExtendedFreeSurferSegmentation: Combined segmentation
    """
    # Compute overlap
    overlap = (seg1.mri.data > 0) & (seg2.mri.data > 0)

    # Offset seg2 labels
    seg2_offset = np.where(seg2.mri.data > 0, seg2.mri.data + offset, 0)

    combined_labels = seg1.mri.data
    combined_labels += seg2_offset

    # Overwrite overlap depending on priority
    if priority == "right":
        combined_labels[overlap] = seg2_offset[overlap]
    elif priority == "left":
        combined_labels[overlap] = seg1.mri.data[overlap]
    else:
        raise Exception

    return ExtendedFreeSurferSegmentation(
        MRIData(combined_labels, seg1.mri.affine),
        seg1.lut,
    )


def segment_dura(segmentation, looklocker, mixed) -> ExtendedFreeSurferSegmentation:
    # Make anatomically based estimate of where the dura should be
    ribbon = make_ribbon(segmentation)

    dura_mask = extract_dura_threshold(ribbon, looklocker, mixed)

    # Extrapolate rois from tissue csf segmentation
    add_dura_to_segmentation(dura_mask, segmentation)

    return segmentation


def make_ribbon(segmentation: ExtendedFreeSurferSegmentation) -> np.ndarray:
    parenchyma_mask = np.where(
        np.logical_and(segmentation.mri.data > 0, segmentation.mri.data < 10000),
        1,
        0,
    )
    parenchyma = binary_erosion(parenchyma_mask, iterations=1)
    # Make ribbon outside CSF spaces
    filled_labels = binary_dilation(
        segmentation.mri.data,
        iterations=6,
    )  # Fill the entire mask
    inner_boundary = binary_erosion(filled_labels, iterations=6).astype(int)
    outer_boundary = binary_dilation(inner_boundary, iterations=3).astype(int)
    outer_ribbon = (outer_boundary - inner_boundary) * (~parenchyma)
    ribbon = outer_ribbon
    ribbon[..., :136] = 0  # Manual adjustment to remove bottom edge of CSF

    labels = measure.label(ribbon, connectivity=1)

    # Only keep biggest region, which should be the dura.
    # This is a bit hacky but it seems to work well in practice
    # Necessary to account for missing regions in the segmentation
    # which can cause the ribbon to be disconnected
    counts = np.bincount(labels.flat)[1:]
    max_label = np.argmax(counts) + 1
    ribbon = ribbon * (labels == max_label)

    return ribbon


def extract_dura_threshold(
    ribbon: np.ndarray,
    LL: np.ndarray,
    Mixed: np.ndarray,
) -> np.ndarray:
    # Ensure all CSF signal is removed:
    LL_masked = LL * ribbon
    dura_mask = (LL_masked < np.percentile(LL[Mixed > 0], 5)) & (LL_masked > 400)

    return dura_mask.astype(int)


def add_dura_to_segmentation(
    mask: np.ndarray,
    segmentation: ExtendedFreeSurferSegmentation,
) -> ExtendedFreeSurferSegmentation:
    # Find nearest neighbour label from CSF segmentation
    label_grid = np.argwhere(segmentation.mri.data > 0)
    mask_grid = np.argwhere(mask > 0)
    interpolator = NearestNDInterpolator(
        label_grid,
        segmentation.mri.data[*label_grid.T],
    )
    seg_values = interpolator(*mask_grid.T)
    # Dura rois are given ids > 20000
    seg_values = np.where(seg_values > 10000, seg_values + 10000, seg_values + 20000)

    segmentation.mri.data[*mask_grid.T] = seg_values
