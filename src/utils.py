import numpy as np


def mask_rois(arr: np.ndarray, rois: np.ndarray, mask_value: float, base_value: float):
    return np.where(
        np.isin(
            arr,
            rois,
        ),
        mask_value,
        base_value,
    )
