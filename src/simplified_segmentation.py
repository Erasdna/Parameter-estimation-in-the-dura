import numpy as np
from src.utils import mask_rois
from mritk.segmentation import (
    ExtendedFreeSurferSegmentation,
    procedural_freesurfer_lut,
    Segmentation,
)
from mritk import MRIData
import pandas as pd


def make_simplified_segmentation(
    segmentation: ExtendedFreeSurferSegmentation,
) -> np.ndarray:
    # Remove ventricles (not relevant for dura):
    reference_segmentation = mask_rois(
        segmentation.mri.data,
        10000 + np.array([4, 5, 14, 15, 43, 44]),
        0,
        segmentation.mri.data,
    )

    # Cerebellum and brainstem
    cerebellum_brainstem = np.array([7, 8, 46, 47, 16], dtype=int)

    # Add cerebellar CSF rois to simplified segmentation
    simplified_segmentation = mask_rois(
        reference_segmentation,
        10000 + cerebellum_brainstem,
        1,
        0,
    )
    # Add cerebellar dura rois to simplified segmentation
    simplified_segmentation = mask_rois(
        reference_segmentation,
        20000 + cerebellum_brainstem,
        10,
        simplified_segmentation,
    )

    # Remove cerebellar rois from csf
    reference_segmentation = mask_rois(
        reference_segmentation,
        10000 + cerebellum_brainstem,
        0,
        reference_segmentation,
    )
    reference_segmentation = mask_rois(
        reference_segmentation,
        20000 + cerebellum_brainstem,
        0,
        reference_segmentation,
    )

    # Split csf and dura into front/back and upper/lower parts based on y and z coordinates
    y_mid = 255
    z__mid = 313

    def assign_split():
        arr = np.zeros_like(reference_segmentation)
        arr[:, :y_mid, :z__mid] = 2  # front lower
        arr[:, :y_mid, z__mid:] = 3  # front upper
        arr[:, y_mid:, :z__mid] = 5  # back lower
        arr[:, y_mid:, z__mid:] = 4  # back upper
        return arr.astype(simplified_segmentation.dtype)

    simplified_segmentation += (
        10 * (reference_segmentation >= 20000) * assign_split()
    ) + (
        (reference_segmentation >= 10000) * (reference_segmentation < 20000)
    ) * assign_split()

    naming = {
        1: "cerebellum_brainstem",
        2: "anterior-inferior",  # front lower
        3: "anterior-superior",  # front upper
        4: "posterior-superior",  # back upper
        5: "posterior-inferior",  # back lower
    }

    lut = procedural_freesurfer_lut(list(naming.keys()), list(naming.values()))
    dura_lut = lut.copy()
    dura_lut.index *= 10
    dura_lut["description"] = dura_lut["description"] + "-dura"
    lut = pd.concat([lut, dura_lut])

    return Segmentation(
        MRIData(simplified_segmentation.astype(float), segmentation.mri.affine),
        pd.DataFrame(lut),
    )
