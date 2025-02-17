#!/usr/bin/env python
"""
Script to perform gradient analysis comparing pre and post FUS.
I think the idea here is to assess whether preoperative alterations (within the sample) may constrain surgery-induced changes.

If preoperative anatomy constrains improvement, you might expect:

Less deviation from baseline (i.e., stable gradients) in poor responders.
Larger gradient shifts in patients with greater flexibility in brain organization.
If anatomy constrains change, the degree of correlation should be higher in those with poorer response?

In order to capture the effects of targeting, we can use the immediate score to dictate accuracy? However, the flexibility of the brain to respond may mediate longer term effcets.

# I have 3 sessions here. What is the best way to approach this?

Option 1: perfrom leave-one-out alignment using baseline data.
    Essentially, construct a baseline (pre-op) template from all but that subject, perform DDT.
    Compute the individuals gradient (perhaps all gradients here) and then align them to the group with Procrustes - calculate the 
    difference between aligned as change (can do this for each hemisphere). Could also compute spin permutations between the maps.
    
"""

# Do not have to perform DDT as this is a problem specific to the biases in short-long range connections of probabilistic tractography.

import os
import glob
import numpy as np
import pandas as pd
from nilearn import plotting, datasets
from brainspace.utils.parcellation import map_to_labels
from brainspace.gradient import GradientMaps
from brainspace.datasets import load_conte69, load_parcellation
from brainspace.plotting import plot_hemispheres
from brainspace.gradient.alignment import procrustes_alignment
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import ttest_ind, ks_2samp, mannwhitneyu, shapiro
from scipy.spatial.distance import pdist, squareform
from statsmodels.stats.multitest import multipletests
from brainspace.datasets import load_group_fc
from brainspace.null_models import SpinPermutations
from scipy.stats import spearmanr

################################
# 1. Load parcellation and groupings
################################

# Load the Schaefer parcellation (400 parcels)
labeling = load_parcellation('schaefer', scale=400, join=True)
mask = labeling != 0
lut_file = "/Users/a9ws/imaging/datasets/mrgfus/derivatives/study_files/parcellations/lut/Schaefer2018_400Parcels_7Networks_order.lut"

# Read subject data (includes subject_id, session, age, sex, clinical scores, etc.)
data_path = "/Users/a9ws/imaging/datasets/mrgfus/derivatives/study_files/gradients/groupings.csv"
data = pd.read_csv(data_path)

# Ensure subject_id is a string
data['subject_id'] = data['subject_id'].astype(str)

# Organize sessions (assuming a 'session' column exists with values like 'ses-01', 'ses-02', 'ses-03')
ses01_ids = data[data['session'] == 'ses-01']['subject_id'].tolist()
ses02_ids = data[data['session'] == 'ses-02']['subject_id'].tolist()
ses03_ids = data[data['session'] == 'ses-03']['subject_id'].tolist()

################################
# 2. Load connectivity matrices
################################

base_dir = "/Users/a9ws/imaging/datasets/mrgfus/derivatives/data/mind_networks"

def load_connectivity_matrix(subj_id, session):
    """
    Load subject-specific connectivity matrix for a given session.
    """
    if not subj_id.startswith("sub-"):
        subj_id = "sub-" + subj_id
    file_path = os.path.join(base_dir, f"{subj_id}_{session}_MIND.csv")
    if os.path.isfile(file_path):
        return np.loadtxt(file_path, delimiter=',')
    else:
        print(f"Connectivity matrix not found for {subj_id}, {session} at {file_path}")
        return None

# Load connectivity matrices for each session
ses01_matrices = {subj: load_connectivity_matrix(subj, "ses-01") for subj in ses01_ids}
ses02_matrices = {subj: load_connectivity_matrix(subj, "ses-02") for subj in ses02_ids}
ses03_matrices = {subj: load_connectivity_matrix(subj, "ses-03") for subj in ses03_ids}

# Remove None entries (subjects with missing matrices)
ses01_matrices = {k: v for k, v in ses01_matrices.items() if v is not None}
ses02_matrices = {k: v for k, v in ses02_matrices.items() if v is not None}
ses03_matrices = {k: v for k, v in ses03_matrices.items() if v is not None}

################################
# 3. Organizing a DataFrame to Track Sessions and Scores
################################

# Define clinical metrics of interest
clinical_vars = ["age", "sex", "group", "UPDRS", "pain_score"]

# Create a long-format DataFrame for easier tracking
session_data = []
for session, matrices in zip(["ses-01", "ses-02", "ses-03"], 
                             [ses01_matrices, ses02_matrices, ses03_matrices]):
    for subj, matrix in matrices.items():
        row = {"subject_id": subj, "session": session, "matrix": matrix}
        for var in clinical_vars:
            row[var] = data.loc[data['subject_id'] == subj, var].values[0]
        session_data.append(row)

# Convert to DataFrame
df_sessions = pd.DataFrame(session_data)

# Computes reference gradients for each patient using leave-one-out construction
# Aligns each patients gradients to their reference gradient using procrustes alignment

from brainspace.gradient import GradientMaps
from brainspace.gradient.alignment import procrustes_alignment
import numpy as np

def compute_group_reference_gradient(exclude_subject, ses01_matrices):
    """
    Computes the reference gradient using all ses01 subjects except `exclude_subject`.
    """
    matrices = [mat for subj, mat in ses01_matrices.items() if subj != exclude_subject]
    
    if len(matrices) == 0:
        raise ValueError(f"No other subjects available to compute reference for {exclude_subject}")
    
    group_avg = np.mean(matrices, axis=0)  # Compute group average connectivity matrix
    gm = GradientMaps(n_components=10, approach='dm', kernel='normalized_angle')  # Example parameters
    gm.fit(group_avg)  # Fit gradients to group reference
    return gm.gradients_

def align_patient_gradients(patient_id, ses01_matrices, ses02_matrices, ses03_matrices, reference_gradient):
    """
    Computes patient-specific gradients for ses01, ses02, ses03 and aligns them to the reference gradient.
    """
    gm = GradientMaps(n_components=10, approach='dm', kernel='normalized_angle')
    
    matrices = {
        "ses01": ses01_matrices.get(patient_id),
        "ses02": ses02_matrices.get(patient_id),
        "ses03": ses03_matrices.get(patient_id)
    }
    
    aligned_gradients = {}
    for session, matrix in matrices.items():
        if matrix is not None:
            gm.fit(matrix)  # Compute patient-specific gradients
            aligned = procrustes_alignment(gm.gradients_, reference_gradient)  # Align to reference
            aligned_gradients[session] = aligned

    return aligned_gradients

# Compute and align for every subject
aligned_gradients_per_subject = {}

for subject in ses01_matrices.keys():
    print(f"Processing subject: {subject}")
    
    # Compute subject-specific reference gradient
    reference_gradient = compute_group_reference_gradient(subject, ses01_matrices)
    
    # Align subject's gradients to their reference
    aligned_gradients_per_subject[subject] = align_patient_gradients(subject, ses01_matrices, ses02_matrices, ses03_matrices, reference_gradient)

print("All subjects processed.")


# Global (lh+rh) correlations of gradients can be performed as well as sorting the primary affected hemisphere:
# We have df column:LESION_SIDE; specifies the hemisphere of the lesion on the lh or rh - need some way to be comparing the gradient for left and right FUS lesions.
# Then also need to compute region and network-wide differences in relation.












# To compare group maps of gradients at different time points

# Convert the connectivity matrices (currently stored in dictionaries) into 3D arrays.
def mean_connectivity_matrices(mat_dict):
    """Stack and average matrices from a dictionary into a 2D mean connectivity matrix."""
    matrices = list(mat_dict.values())
    if matrices:
        stacked_matrices = np.stack(matrices, axis=2)  # Shape (n, n, s)
        return np.mean(stacked_matrices, axis=2)  # Average over subjects (s)
    else:
        raise ValueError("No matrices found in the provided dictionary.")

ses01_group = mean_connectivity_matrices(ses01_matrices)
ses02_group = mean_connectivity_matrices(ses02_matrices)
ses03_group = mean_connectivity_matrices(ses03_matrices)

################################
# 5. Map the session PG's
################################

# Load hemisphere labels from the Schaefer parcellation LUT file
lut_df = pd.read_csv(lut_file, sep="\s+", header=None, names=["Index", "Label"])
parcel_labels = lut_df["Label"].astype(str).tolist()

# Function to compute gradients
def compute_gradients(matrix, n_components=10):
    """
    Compute diffusion map gradients for a given connectivity matrix.
    Returns the GradientMaps object.
    """
    gm = GradientMaps(approach='dm', n_components=10, random_state=0)
    gm.fit(matrix)
    return gm  # Return the full object, not just the gradients array

# Compute gradients for each session
ses01_grad = compute_gradients(ses01_group, n_components=10)
ses01_grad = ses01_grad.gradients_
ses01_pg = map_to_labels(ses01_grad[:, 0], labeling, mask=mask, fill=np.nan)

ses02_grad = compute_gradients(ses02_group, n_components=10)
ses02_grad = ses02_grad.gradients_
ses02_pg = map_to_labels(ses02_grad[:, 0], labeling, mask=mask, fill=np.nan)

ses03_grad = compute_gradients(ses03_group, n_components=10)
ses03_grad = ses01_grad.gradients_
ses03_pg = map_to_labels(ses03_grad[:, 0], labeling, mask=mask, fill=np.nan)

# Plot the principal gradients
surf_lh, surf_rh = load_conte69()
plot_hemispheres(
    surf_lh, surf_rh,
    array_name=[ses01_pg, ses02_pg, ses03_pg],
    size=(1200, 600),
    cmap='RdBu',
    color_bar=True,
    label_text=['Baseline G1', '6-month G1', '12-month G1'],
    zoom=1.05
)
plt.show()

# Spin permutations

from brainspace.null_models import SpinPermutations
from brainspace.datasets import load_conte69
from scipy.stats import spearmanr
# Load sphere surfaces for spin testing
sphere_lh, sphere_rh = load_conte69(as_sphere=True)

# Ensure the shapes of the sphere surfaces match the gradient mappings
n_vertices_lh = sphere_lh.points.shape[0]
n_vertices_rh = sphere_rh.points.shape[0]

# Ensure gradient mappings match the surface resolution
assert len(ses01_pg) == n_vertices_lh + n_vertices_rh, "Ses01 G1 mapping size mismatch."
assert len(ses02_pg) == n_vertices_lh + n_vertices_rh, "Ses02 G1 mapping size mismatch."
assert len(ses03_pg) == n_vertices_lh + n_vertices_rh, "Ses03 G1 mapping size mismatch."

# Spin Permutations
n_rand = 5000  # Number of random permutations
sp = SpinPermutations(n_rep=n_rand, random_state=0)
sp.fit(sphere_lh, points_rh=sphere_rh)  # Fit to left and right hemisphere spheres

# Rotate ses01 gradient mappings
ses01_rotated = np.hstack(
    sp.randomize(ses01_pg[:n_vertices_lh], ses01_pg[n_vertices_lh:])
)

# Compute observed correlation between gradients
mask = ~np.isnan(ses03_pg) & ~np.isnan(ses01_pg)
r_obs, pv_obs = spearmanr(ses03_pg[mask], ses01_pg[mask])

# Compute correlations for permuted gradients
r_spin = np.empty(n_rand)
for i, perm in enumerate(ses01_rotated):
    mask_rot = mask & ~np.isnan(perm)
    r_spin[i] = spearmanr(ses03_pg[mask_rot], perm[mask_rot])[0]

# Compute spin-test p-value
pv_spin = np.mean(np.abs(r_spin) >= np.abs(r_obs))

# Plot Null Distribution
fig, ax = plt.subplots(figsize=(7, 5))
ax.hist(r_spin, bins=25, density=True, alpha=0.5, color=(.8, .8, .8), label='Null Distribution')
ax.axvline(r_obs, lw=2, ls='--', color='k', label=f'Observed r = {r_obs:.3f}')
ax.set_title('Spin Permutation Test')
ax.set_xlabel('Spearman Correlation')
ax.set_ylabel('Density')
ax.legend()

plt.show()

# Print Results
print(f"Observed Correlation (r): {r_obs:.3f}")
print(f"Observed P-Value: {pv_obs:.5e}")
print(f"Spin Permutation P-Value: {pv_spin:.5e}")

# Plot the first combined gradient
plot_hemispheres(
    surf_lh, surf_rh,
    array_name=[ses01_pg, ses03_pg],
    size=(1200, 600),
    cmap='RdBu',
    color_bar=True,
    label_text=['Ses01 PG', 'Ses03 PG'],
    zoom=1.05
)
plt.show()

