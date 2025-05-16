# File: segment_nifti_pca.py

import os
import numpy as np
import nibabel as nib
from sklearn.decomposition import PCA

def segment_nifti_pca(input_nifti_path, output_dir, N):
    # Load NIfTI mask
    nii = nib.load(input_nifti_path)
    vol = nii.get_fdata()
    affine = nii.affine
    header = nii.header

    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Get coordinates of non-zero voxels
    coords = np.array(np.nonzero(vol)).T  # shape: (num_voxels, 3)

    # Perform PCA
    pca = PCA(n_components=3)
    projected = pca.fit_transform(coords)  # shape: (num_voxels, 3)

    for pc in range(3):  # PC1, PC2, PC3
        proj_vals = projected[:, pc]
        edges = np.linspace(proj_vals.min(), proj_vals.max(), N + 1)

        for i in range(N):
            in_bin = (proj_vals >= edges[i]) & (proj_vals < edges[i + 1])
            selected_coords = coords[in_bin]

            # Create binary mask
            mask = np.zeros(vol.shape, dtype=np.uint8)
            mask[tuple(selected_coords.T)] = 1

            # Save to NIfTI
            out_nii = nib.Nifti1Image(mask, affine, header)
            out_path = os.path.join(output_dir, f'segment_PC{pc+1}_{i+1:02d}.nii.gz')
            nib.save(out_nii, out_path)

    print(f"Saved {N} segmented masks along each of 3 PCA axes to {output_dir}")
