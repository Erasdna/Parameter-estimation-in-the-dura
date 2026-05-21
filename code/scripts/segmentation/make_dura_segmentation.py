from snakemake.script import snakemake
import nibabel as nib
from pathlib import Path
from mritk.segmentation import (
    FreeSurferSegmentation,
    ExtendedFreeSurferSegmentation,
)
from mritk import MRIData
import numpy as np
from scipy.ndimage import binary_dilation
from src.dura_segmentation import combine_segmentations, segment_dura
from src.simplified_segmentation import make_simplified_segmentation

wmparc_path = Path(snakemake.input.wmparc)
csf_path = Path(snakemake.input.csf)

wmparc = FreeSurferSegmentation.from_file(wmparc_path)
csf = FreeSurferSegmentation.from_file(csf_path)

# Combine parenchymal an CSF segmentations
combined = combine_segmentations(wmparc, csf, offset=10000, priority="right")

# Run dura segmentation, using the mixed and look-locker sequences for thresholding
mixed = nib.load(snakemake.input.mixed).get_fdata().astype(int)
looklocker = nib.load(snakemake.input.looklocker).get_fdata().astype(int)

dura_segmentation = segment_dura(combined, looklocker, mixed)
dura_segmentation.save(
    Path(snakemake.params.dura_segmentation),
)

# Make simplified segmentation, splitting the brain into 5 main regions
simplified_segmentation = make_simplified_segmentation(dura_segmentation)
simplified_segmentation.save(
    Path(snakemake.params.simplified_segmentation),
)

# Extract outer CSF layer
dura_mask = np.where(dura_segmentation.mri.data >= 20000, 1, 0)
dilated_dura = binary_dilation(dura_mask, iterations=3)
simplified_segmentation_layer_labels = simplified_segmentation.mri.data * dilated_dura
simplified_segmentation_layer = ExtendedFreeSurferSegmentation(
    MRIData(simplified_segmentation_layer_labels, simplified_segmentation.mri.affine),
    simplified_segmentation.lut,
)
simplified_segmentation_layer.save(
    Path(snakemake.params.simplified_segmentation_layer),
)
