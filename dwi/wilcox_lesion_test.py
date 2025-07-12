# Python code to implement Wilcoxon voxel testing for MRgFuS lesions

import numpy as np
import nibabel as nib
import os
from glob import glob
from scipy.stats import wilcoxon
from tqdm import tqdm

# === Paths ===
lesion_dir = '~/imaging/mrgfus/derivatives/lesions/data'
score_file = '~/imaging/mrgfus/derivatives/study_files/imm_scores.txt'

# === Load Clinical Scores ===
clinical_scores = np.loadtxt(score_file)
group_mean = np.mean(clinical_scores)
n_subjects = len(clinical_scores)

# === Load Lesion Masks ===
lesion_paths = sorted(glob(os.path.join(lesion_dir, 'sub-*', '*_lesion.nii.gz')))
assert len(lesion_paths) == n_subjects, "Mismatch between lesion files and clinical scores"

# Load all lesions into a 4D array
print("Loading lesion masks...")
lesion_data = []
for path in lesion_paths:
    img = nib.load(path)
    data = img.get_fdata().astype(bool)
    lesion_data.append(data)

lesion_data = np.stack(lesion_data, axis=-1)  # Shape: (X, Y, Z, N)
ref_img = nib.load(lesion_paths[0])  # For header/affine

# === Compute Voxel-Wise Coverage ===
print("Computing voxel coverage...")
coverage = np.sum(lesion_data, axis=-1)  # Shape: (X, Y, Z)
threshold = int(0.2 * n_subjects)
included_voxels = coverage >= threshold

# === Run Voxel-Wise Wilcoxon Tests ===
print("Running Wilcoxon tests...")
z_map = np.full(lesion_data.shape[:3], np.nan)
p_map = np.full(lesion_data.shape[:3], np.nan)

voxel_indices = np.argwhere(included_voxels)

for idx in tqdm(voxel_indices, desc="Processing voxels"):
    x, y, z = idx
    voxel_mask = lesion_data[x, y, z, :]
    subj_scores = clinical_scores[voxel_mask]
    
    if len(subj_scores) < 2:
        continue  # Not enough data for test
    
    try:
        stat, p = wilcoxon(subj_scores - group_mean, alternative='greater')
        z_map[x, y, z] = stat
        p_map[x, y, z] = p
    except:
        continue  # Handle rare exceptions gracefully

# === Save Outputs ===
print("Saving results...")
z_img = nib.Nifti1Image(z_map, ref_img.affine, ref_img.header)
p_img = nib.Nifti1Image(p_map, ref_img.affine, ref_img.header)

nib.save(z_img, '~/imaging/mrgfus/derivatives/lesions/sweetspot_zmap.nii.gz')
nib.save(p_img, '~/imaging/mrgfus/derivatives/lesions/sweetspot_pmap.nii.gz')

from statsmodels.stats.multitest import multipletests

# === Flatten and Extract p-values for included voxels only ===
print("Applying FDR correction...")
p_values = p_map[included_voxels]

# FDR correction (Benjamini-Hochberg)
reject_fdr, pvals_corrected_fdr, _, _ = multipletests(p_values, alpha=0.05, method='fdr_bh')

# Create maps to hold corrected results
fdr_corrected_p_map = np.full(p_map.shape, np.nan)
fdr_reject_map = np.zeros(p_map.shape, dtype=bool)

# Reinsert corrected p-values and reject mask into voxel-wise maps
fdr_corrected_p_map[included_voxels] = pvals_corrected_fdr
fdr_reject_map[included_voxels] = reject_fdr

# Output FDR-corrected p-value map
fdr_p_img = nib.Nifti1Image(fdr_corrected_p_map, ref_img.affine, ref_img.header)
nib.save(fdr_p_img, '~/imaging/mrgfus/derivatives/lesions/sweetspot_pmap_FDR_corrected.nii.gz')

# Output binary sweetspot mask (thresholded at FDR < 0.05)
fdr_binary_img = nib.Nifti1Image(fdr_reject_map.astype(np.uint8), ref_img.affine, ref_img.header)
nib.save(fdr_binary_img, '~/imaging/mrgfus/derivatives/lesions/sweetspot_mask_FDR_p05.nii.gz')
