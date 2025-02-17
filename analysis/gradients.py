#!/usr/bin/env python
"""
Script to perform gradient analysis comparing pre and post FUS.
I think the idea here is to assess whether preoperative alterations (within the sample) may constrain surgery-induced changes.

"""

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

# Read groupings information (this CSV contains subject_id, site, age, sex, group, etc.)
data_path = "/Users/a9ws/imaging/datasets/pd/derivatives/study_files/gradients/groupings.csv"
data = pd.read_csv(data_path)
# Ensure subject_id is string and set it as index for later lookup.
data['subject_id'] = data['subject_id'].astype(str)
data.set_index("subject_id", inplace=True)

# Filter groups:
# For example, group 1 = patients (60d) and group 3 = healthy controls.
patients_ids = data.loc[data['group'] == 1].index.values  # Patients
older_hc_ids = data.loc[data['group'] == 3].index.values    # Standard healthy controls
young_hc_ids = data.loc[data['group'] == 4].index.values      # Young healthy controls



# Implement sorting (need some way to be comparing the gradient for left and right FUS lesions)


################################
# 1.5. Load ROI centroids for distance-dependent thresholding
################################

centroids_file = "/Users/a9ws/imaging/datasets/pd/derivatives/study_files/parcellations/centroids/Schaefer2018_400Parcels_7Networks_order_FSLMNI152_1mm.Centroid_RAS.csv"
centroids_df = pd.read_csv(centroids_file, sep=",")
roi_centers = centroids_df[['R', 'A', 'S']].to_numpy()  # shape: (400, 3)

################################
# 2. Load connectivity matrices
################################

base_dir = "/Users/a9ws/imaging/datasets/mrgfus/derivatives/data/mind_networks"

# Need to read in the sessions seperately...

def load_connectivity_matrix(subj_id, filename="Schaefer_400_tck_gradmat.csv"):
    """
    Given a subject id, construct the expected file path and load the connectivity matrix.
    """
    if not subj_id.startswith("sub-"):
        subj_id = "sub-" + subj_id
    file_path = os.path.join(base_dir, filename)
    if os.path.isfile(file_path):
        return np.loadtxt(file_path, delimiter=',')
    else:
        print(f"Connectivity matrix not found for {subj_id} at {file_path}")
        return None

# Load connectivity matrices for each timepoint
ses_01_matrices = {}
for subj in patients_ids:
    mat = load_connectivity_matrix(subj, filename="${sub}_ses-01_MIND.csv")
    if mat is not None:
        patients_matrices[subj] = mat

ses_02_matrices = {}
for subj in older_hc_ids:
    mat = load_connectivity_matrix(subj, filename="${sub}_ses-02_MIND.csv")
    if mat is not None:
        older_hc_matrices[subj] = mat

ses_03_matrices = {}
for subj in young_hc_ids:
    mat = load_connectivity_matrix(subj, filename="${sub}_ses-03_MIND.csv")
    if mat is not None:
        young_hc_matrices[subj] = mat

################################
# 3. Define the distance-dependent thresholding function
################################

def distance_thresholding(A, dist, hemiid, nbins):
    '''
    Input:
    - A: n x n x s 
        connection strengths between n edges for s subjects.
    - dist: n x n 
        edge distances
    - hemiid: n x 1 
        vector indicating the hemisphere of each region (e.g., 0 for left, 1 for right)
    - nbins (int) 
        number of bins for distance thresholding
    Output:
    - G: n x n 
        group connectivity matrix using distance thresholding
    - Gc: n x n 
        group connectivity matrix using classic consensus thresholding
    '''
    
    n, _, nsub = np.shape(A)    # number nodes (n) and number of subjects (nsub)
    C = np.sum(A > 0, axis=2)     # consistency (number of subjects with a connection)
    W = np.sum(A, axis=2) / C     # average weight
    W[np.isnan(W)] = 0           # remove NaNs
    Grp = np.zeros((n, n, 2))     # will store inter- and intra-hemispheric connections separately
    Gc = Grp.copy()
    
    # Create bins based on edge distances.
    # (Note: np.nonzero(dist) returns indices; you might want to use dist[dist > 0] instead.)
    distbins = np.linspace(np.min(np.nonzero(dist)), np.max(np.nonzero(dist)), nbins + 1)
    distbins[-1] = distbins[-1] + 1
    
    for j in range(2):  # j==0: inter-hemispheric, j==1: intra-hemispheric
        # Create mask for inter- or intra-hemispheric edges.
        if j == 0:  # inter-hemispheric
            d = (hemiid == 0) * (hemiid.T == 1)
            d = np.logical_or(d, d.T)
        else:       # intra-hemispheric
            d = ((hemiid == 0) * (hemiid.T == 0)) | ((hemiid == 1) * (hemiid.T == 1))
            d = np.logical_or(d, d.T)
        
        # Compute a weighted distance matrix for each subject
        D = np.empty_like(A)
        upper_mask = dist * np.triu(d)
        for i in range(nsub):
            binary_A = A > 0
            D[:, :, i] = binary_A[:, :, i] * upper_mask
        D = D[np.nonzero(D)].flatten()
        
        # Determine target number of edges to retain
        tgt = len(D) / nsub
        G_flat = np.zeros((n*n))
        
        # Process each distance bin
        for ibin in range(nbins):
            # Determine fraction of connections to keep in this bin
            D_lower = D >= distbins[ibin]
            D_upper = D < distbins[ibin + 1]
            frac = round((tgt * np.sum(np.logical_and(D_lower, D_upper))) / len(D))
            
            # Mask out connections not in current bin
            c = np.triu(np.logical_and(dist >= distbins[ibin], dist < distbins[ibin + 1])) * C * d
            
            # Rank connections (descending order) based on consistency
            sorted_idx = np.argsort(c, axis=None)[::-1]
            
            # Keep connections until reaching the target fraction for this bin
            G_flat[sorted_idx[:frac]] = 1
          
        # Reshape and store the distance-thresholded connectivity for this hemisphere grouping
        Grp[:, :, j] = np.reshape(G_flat, [n, n])
        
        # For classic consensus thresholding: select top connections based on average weight
        w = W * np.triu(d)
        idx = np.argsort(w, axis=None)[::-1]
        w_flat = np.zeros((n*n))
        w_flat[idx[:np.count_nonzero(G_flat)]] = 1
        Gc[:, :, j] = np.reshape(w_flat, [n, n])
    
    # Combine inter- and intra-hemispheric matrices and symmetrize the final graph
    G = np.sum(Grp, axis=2)
    G = np.maximum(G, G.transpose())
    
    Gc = np.sum(Gc, axis=2)
    Gc = np.maximum(Gc, Gc.transpose())
        
    return G, Gc

################################
# 4. Perform distance-dependent thresholding for each group
################################

# I have 3 sessions here. What is the best way to approach this?

Option 1: perfrom leave-one-out alignment using baseline data.
    Essentially, construct a baseline (pre-op) template from all but that subject, perform DDT.
    Compute the individuals gradient (perhaps all gradients here) and then align them to the group with Procrustes - calculate the 
    difference between aligned as change (can do this for each hemisphere). Could also compute spin permutations between the maps.
    If anatomy constrains change, the degree of correlation should be higher in those with poorer response?




# Compute the pairwise Euclidean distance between ROI centroids
dist = squareform(pdist(roi_centers, metric='euclidean'))

# Define hemisphere IDs.
# Here, we assume the first half of parcels are left hemisphere (0) and the second half are right (1).
n_rois = roi_centers.shape[0]
hemiid = np.zeros((n_rois, 1), dtype=int)
hemiid[n_rois // 2:] = 1  # first half = 0 (left), second half = 1 (right)

# Define number of distance bins (adjust as needed)
nbins = 5

# Convert the connectivity matrices (currently stored in dictionaries) into 3D arrays.
def stack_connectivity_matrices(mat_dict):
    """Stack matrices from a dictionary into a 3D array (n x n x s)."""
    matrices = list(mat_dict.values())
    if matrices:
        return np.stack(matrices, axis=2)
    else:
        raise ValueError("No matrices found in the provided dictionary.")

patients_A = stack_connectivity_matrices(patients_matrices)
older_hc_A = stack_connectivity_matrices(older_hc_matrices)
young_hc_A = stack_connectivity_matrices(young_hc_matrices)

# Apply distance-dependent thresholding for each group (_G matrices are binary masks)
patients_G, patients_Gc = distance_thresholding(patients_A, dist, hemiid, nbins)
older_hc_G, older_hc_Gc = distance_thresholding(older_hc_A, dist, hemiid, nbins)
young_hc_G, young_hc_Gc = distance_thresholding(young_hc_A, dist, hemiid, nbins)

# Extract the mean streamlines across the groups
young_hc_mean_streamlines = np.mean(young_hc_A, axis=2)
# Apply the mask
thresholded_young_hc_mean_streamlines = young_hc_mean_streamlines * young_hc_G

older_hc_mean_streamlines = np.mean(older_hc_A, axis=2)
thresholded_older_hc_mean_streamlines = older_hc_mean_streamlines * older_hc_G

patient_mean_streamlines = np.mean(patients_A, axis=2)
thresholded_patient_mean_streamlines = patient_mean_streamlines * patients_G

################################
# 5. Compute the healthy young template and gradients
################################

# I want to compute the healthy young template here but want the left and right hemisphere gradients to be computed seperately and then
# combined afterwards to create gradients.

def split_hemispheres(matrix, labels):
    """
    Splits a connectivity matrix into left and right hemisphere matrices.
    
    Parameters:
    - matrix: (n_rois, n_rois) connectivity matrix
    - labels: (n_rois,) list of parcel labels indicating hemisphere (e.g., 'LH', 'RH')
    
    Returns:
    - matrix_LH: (n_LH, n_LH) left hemisphere connectivity matrix
    - matrix_RH: (n_RH, n_RH) right hemisphere connectivity matrix
    """
    lh_indices = [i for i, label in enumerate(labels) if 'LH' in label]
    rh_indices = [i for i, label in enumerate(labels) if 'RH' in label]

    matrix_LH = matrix[np.ix_(lh_indices, lh_indices)]
    matrix_RH = matrix[np.ix_(rh_indices, rh_indices)]
    
    return matrix_LH, matrix_RH, lh_indices, rh_indices

# Load hemisphere labels from the Schaefer parcellation LUT file
lut_df = pd.read_csv(lut_file, sep="\s+", header=None, names=["Index", "Label"])
parcel_labels = lut_df["Label"].astype(str).tolist()

# Split the young healthy template into LH and RH
young_hc_template_LH, young_hc_template_RH, lh_indices, rh_indices = split_hemispheres(thresholded_young_hc_mean_streamlines, parcel_labels)

# Function to compute gradients (returns the full GradientMaps object)
def compute_gradients(matrix, n_components=3):
    """
    Compute diffusion map gradients for a given connectivity matrix.
    Returns the GradientMaps object.
    """
    gm = GradientMaps(approach='dm', n_components=n_components, random_state=0)
    gm.fit(matrix)
    return gm  # Return the full object, not just the gradients array

# Compute gradients separately for LH and RH
young_hc_gm_LH = compute_gradients(young_hc_template_LH, n_components=3)
young_hc_gm_RH = compute_gradients(young_hc_template_RH, n_components=3)

# Extract gradients from the objects
young_hc_gradients_LH = young_hc_gm_LH.gradients_  # Shape: (n_LH_parcels, 3)
young_hc_gradients_RH = young_hc_gm_RH.gradients_  # Shape: (n_RH_parcels, 3)

# Initialize a full 400-parcel array for combined gradients
young_hc_gradients_combined = np.zeros((len(parcel_labels), 3))  # Shape: (400, 3)

# Assign LH and RH gradients into their respective positions
young_hc_gradients_combined[lh_indices, :] = young_hc_gradients_LH
young_hc_gradients_combined[rh_indices, :] = young_hc_gradients_RH

# Load Conte69 surface
surf_lh, surf_rh = load_conte69()

# Map the three gradients onto the brain surface
young_hc_g1 = map_to_labels(young_hc_gradients_combined[:, 0], labeling, mask=mask, fill=np.nan)
young_hc_g2 = map_to_labels(young_hc_gradients_combined[:, 1], labeling, mask=mask, fill=np.nan)
young_hc_g3 = map_to_labels(young_hc_gradients_combined[:, 2], labeling, mask=mask, fill=np.nan)

# Plot the first combined gradient
plot_hemispheres(
    surf_lh, surf_rh,
    array_name=[young_hc_g1, young_hc_g2, young_hc_g3],
    size=(1200, 600),
    cmap='RdBu',
    color_bar=True,
    label_text=['Primary Gradient', 'Secondary Gradient', 'Tertiary Gradient'],
    zoom=1.05
)
plt.show()

# HERE IS WHERE IMPLEMENTATION MAY CHANGE
# COMPUTE INDIVIDUAL GRADIENTS (LEFT & RIGHT SEPERATELY AND THEN COMBINE) AND ALIGN TO YOUNG HEALTHY




################################
# 6. Compute connectome perturbations
################################

def compute_perturbed_gradients(subject_matrix, young_hc_template, lh_indices, rh_indices, n_components=3):
    """
    Computes perturbed gradients by averaging the healthy young template (distance-thresholded)
    with the subject's connectivity matrix and recalculating gradients.

    Parameters:
    - subject_matrix: (400, 400) Subject's processed connectivity matrix.
    - young_hc_template: (400, 400) Group representative connectivity matrix (e.g., thresholded_young_hc_mean_streamlines).
    - lh_indices, rh_indices: Hemisphere indices for splitting and recombining gradients.
    - n_components: Number of gradient components to compute.

    Returns:
    - perturbed_gradients_combined: (400, n_components) Combined gradient components after perturbation.
    """
    # Create new mean connectome by averaging the healthy template and the subject matrix.
    new_mean_connectome = np.mean([young_hc_template, subject_matrix], axis=0)

    # Split into left and right hemispheres
    new_mean_LH, new_mean_RH, _, _ = split_hemispheres(new_mean_connectome, parcel_labels)

    # Compute gradients separately for LH and RH
    perturbed_gm_LH = compute_gradients(new_mean_LH, n_components=n_components)
    
    perturbed_gm_RH = compute_gradients(new_mean_RH, n_components=n_components)

    # For instance, if you want the first 3 gradients:
    perturbed_gradients_LH = perturbed_gm_LH.gradients_[:, :3]
    perturbed_gradients_RH = perturbed_gm_RH.gradients_[:, :3]

    # Combine into full-brain gradient matrix
    perturbed_gradients_combined = np.zeros((len(parcel_labels), 3))
    perturbed_gradients_combined[lh_indices, :] = perturbed_gradients_LH
    perturbed_gradients_combined[rh_indices, :] = perturbed_gradients_RH

    return perturbed_gradients_combined

def compute_deviation_index(perturbed_gradients, young_hc_gradients):
    """
    Compute deviation index (DI) using signed log1p transform.

    Parameters:
    - perturbed_gradients: (400, 3) Perturbed gradients after adding a subject
    - young_hc_gradients: (400, 3) Baseline gradients from young HC

    Returns:
    - deviation_index: (400, 3) Signed log1p transformed deviation index
    """
    diff = perturbed_gradients - young_hc_gradients
    deviation_index = np.sign(diff) * np.log1p(np.abs(diff))
    return deviation_index

# Compute DI for older HC's
older_hc_differences = []

for subj in older_hc_matrices:
    perturbed_gradients = compute_perturbed_gradients(older_hc_matrices[subj],
                                                      thresholded_young_hc_mean_streamlines, 
                                                      lh_indices, rh_indices)
    di = compute_deviation_index(perturbed_gradients, young_hc_gradients_combined)
    older_hc_differences.append(di)

# Convert to NumPy array (subjects, 400 parcels, 3 gradients)
older_hc_differences = np.array(older_hc_differences)

# Compute mean and standard deviation across older HC subjects
older_hc_mean = np.mean(older_hc_differences, axis=0)
older_hc_std = np.std(older_hc_differences, axis=0)

# Z-score patient ROIs based on mean/stdev from older HC's
# Store standardized deviation indices for patients
patient_differences_zscored = {}

for subj in patients_matrices:
    perturbed_gradients = compute_perturbed_gradients(patients_matrices[subj], thresholded_young_hc_mean_streamlines, lh_indices, rh_indices)
    di = compute_deviation_index(perturbed_gradients, young_hc_gradients_combined)

    # Standardize using older HC statistics
    z_scored_di = (di - older_hc_mean) / (older_hc_std)
    patient_differences_zscored[subj] = z_scored_di

# Convert dictionary to DataFrame
patient_zscored_df = pd.DataFrame.from_dict(
    {subj: patient_differences_zscored[subj].flatten() for subj in patients_matrices}, 
    orient='index'
)

# Assign meaningful column names
gradient_labels = [f"Gradient{g+1}_Parcel{p+1}" for g in range(3) for p in range(400)]
patient_zscored_df.columns = gradient_labels

# Compute patient ROI eccentricity values based on the older HC manifold origin
# Compute the group manifold origin (mean of the first 3 gradients across older HC)
group_manifold_origin = np.mean(older_hc_differences, axis=0)  # Shape: (400, 3)

def compute_eccentricity(subject_gradients, reference_origin):
    """
    Compute the manifold eccentricity for a patient.

    Parameters:
    - subject_gradients: (400, 3) Subject's deviation indices
    - reference_origin: (400, 3) Group manifold origin (older HC mean)

    Returns:
    - eccentricity: (400,) Euclidean distance per region
    """
    distances = np.linalg.norm(subject_gradients - reference_origin, axis=1)
    return distances  # Shape: (400,)

# Store eccentricity values for all patients
patient_eccentricity = {}

for subj in patients_matrices:
    eccentricity = compute_eccentricity(patient_differences_zscored[subj], group_manifold_origin)
    patient_eccentricity[subj] = eccentricity

# Convert dictionary to DataFrame
eccentricity_df = pd.DataFrame.from_dict(
    {subj: patient_eccentricity[subj] for subj in patients_matrices}, 
    orient='index'
)

# Assign column names for each ROI
eccentricity_df.columns = [f"ROI_{i+1}" for i in range(400)]

# Plot patient averaged eccentricity on the surface
avg_eccentricity = eccentricity_df.mean(axis=0)
std_eccentricity = eccentricity_df.std(axis=0)  # This gives a Series of length equal to the number of ROIs (400)

# Assuming your parcellation 'labeling' and mask are already defined:
avg_eccentricity_mapped = map_to_labels(avg_eccentricity.values, labeling, mask=mask, fill=np.nan)
std_eccentricity_mapped = map_to_labels(std_eccentricity.values, labeling, mask=mask, fill=np.nan)

# Now, plot the values on the brain surfaces using Brainspace's plotting function:
plot_hemispheres(
    surf_lh, surf_rh,
    array_name=[avg_eccentricity_mapped, std_eccentricity_mapped],
    size=(1200, 600),
    cmap='RdBu',
    color_bar=True,
    label_text=['Mean Eccentricity', 'Std Eccentricity'],
    zoom=1.05
)
plt.show()

# 3D scatter plot of average eccentricity values
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.colors as mcolors
import numpy as np

# Suppose you have:
# - coords_3d: (n_rois, 3) array of ROI coordinates (e.g., from group_manifold_origin[:, :3])
# - mean_ecc: (n_rois,) array of mean eccentricity values for each ROI

fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

# Set up colormap normalization based on mean_ecc values
norm = mcolors.Normalize(vmin=avg_eccentricity.min(), vmax=avg_eccentricity.max())
cmap = plt.cm.coolwarm

# Plot each ROI as a point, color-coded by mean eccentricity
scatter = ax.scatter(coords_3d[:, 0], coords_3d[:, 1], coords_3d[:, 2],
                     c=avg_eccentricity, cmap=cmap, norm=norm, s=50, edgecolor='k')

# Define a single "origin" point (using the centroid of all ROI coordinates)
origin = np.mean(coords_3d, axis=0)

# Draw lines from the origin to each ROI, colored by the ROI's eccentricity
for i in range(coords_3d.shape[0]):
    line_color = cmap(norm(avg_eccentricity[i]))
    ax.plot([origin[0], coords_3d[i, 0]],
            [origin[1], coords_3d[i, 1]],
            [origin[2], coords_3d[i, 2]],
            color=line_color, alpha=0.5, linewidth=1)

# Mark the origin with a larger distinct marker
ax.scatter([origin[0]], [origin[1]], [origin[2]], color='black', marker='o', s=150, label='Origin')

# Set the view to bird's-eye (top-down)
ax.view_init(elev=30, azim=-45)

# Instead of removing grid lines and labels, reduce their visual impact:
ax.grid(True, linestyle='--', linewidth=0.5, alpha=0.5)

# Reduce the tick label font size and number of ticks:
ax.tick_params(axis='both', which='major', labelsize=8)

# Optionally, you can set a maximum number of ticks using a locator (example for the x-axis)
from matplotlib.ticker import MaxNLocator
ax.xaxis.set_major_locator(MaxNLocator(nbins=4))
ax.yaxis.set_major_locator(MaxNLocator(nbins=4))
ax.zaxis.set_major_locator(MaxNLocator(nbins=4))

# Keep axis labels but reduce their font size as well
ax.set_xlabel("Gradient 1", fontsize=10)
ax.set_ylabel("Gradient 2", fontsize=10)
ax.set_zlabel("Gradient 3", fontsize=10)
ax.set_title("Mean Eccentricity in 3D Gradient Space (Bird's-eye view)", fontsize=12)

# Add a colorbar
cbar = plt.colorbar(scatter, shrink=0.6, pad=0.1)
cbar.set_label("Mean Eccentricity", fontsize=10)
cbar.ax.tick_params(labelsize=8)

ax.legend(fontsize=10)
plt.show()

# Statistical testing

# Idea: compare the 3 gradients between the older hc and the patient groups.
# Step 1: perform distance dependent thresholding to compute a group matrix from all healthy older controls and all patients.
# Step 2: extract the 3 primary gradients (seperately for the left and right, and then combine)
# Step 3: extract the 3 primary gradients (left and right seperately and then combine for each healthy older control and each patient
# Step 4: use Procrustes alignment to align the individually computed healthy older controls and patients to the group
# Step 5: perform analyses 

# If we want to do this, I believe we should be using procrustes alignment to align each individual matrix to the group one?
# So i think we need to perform distance-dep thresholding on older and PD combined, and then align each individual via Procrustes...

################################
# Compare the 3 gradients between the older HC and the patient groups.
# Step 1: Compute a group connectivity matrix from all older healthy controls and patients.
################################

# Merge older HC and patients dictionaries
combined_matrices = {**older_hc_matrices, **patients_matrices}
combined_A = stack_connectivity_matrices(combined_matrices)
combined_mean_streamlines = np.mean(combined_A, axis=2)

# Apply distance-dependent thresholding on the combined group
combined_G, _ = distance_thresholding(combined_A, dist, hemiid, nbins)
thresholded_combined_mean_streamlines = combined_mean_streamlines * combined_G

################################
# Step 2: Extract the 3 primary gradients for the group template
################################

# Split the combined (group) matrix into left and right hemispheres
group_template_LH, group_template_RH, lh_indices, rh_indices = split_hemispheres(
    thresholded_combined_mean_streamlines, parcel_labels
)

# Compute gradients separately for each hemisphere (we use 3 components)
group_gm_LH = compute_gradients(group_template_LH, n_components=3)
group_gm_RH = compute_gradients(group_template_RH, n_components=3)

# Combine the left and right gradients into a full-brain (400 parcels x 3 gradients) matrix
group_gradients_combined = np.zeros((len(parcel_labels), 3))
group_gradients_combined[lh_indices, :] = group_gm_LH.gradients_[:, :3]
group_gradients_combined[rh_indices, :] = group_gm_RH.gradients_[:, :3]

################################
# Step 3: Define a function to compute subject gradients
################################

def compute_subject_gradients(matrix, labels, lh_idx, rh_idx, n_components=3):
    """
    Compute the subject's primary gradients by splitting the connectivity matrix 
    into hemispheres, computing gradients separately, and then combining them.
    
    Parameters:
      - matrix: (400, 400) subject's connectivity matrix.
      - labels: list of parcel labels (length 400) indicating hemisphere info.
      - lh_idx, rh_idx: indices for left and right hemispheres.
      - n_components: number of gradient components to extract.
      
    Returns:
      - subject_gradients: (400, n_components) gradient matrix.
    """
    # (Optionally, you could also threshold the subject matrix here.)
    subj_LH, subj_RH, _, _ = split_hemispheres(matrix, labels)
    gm_LH = compute_gradients(subj_LH, n_components=n_components)
    gm_RH = compute_gradients(subj_RH, n_components=n_components)
    
    subject_gradients = np.zeros((len(labels), n_components))
    subject_gradients[lh_idx, :] = gm_LH.gradients_[:, :n_components]
    subject_gradients[rh_idx, :] = gm_RH.gradients_[:, :n_components]
    return subject_gradients

################################
# Step 4: Align individual gradients to the group template using Procrustes alignment
################################

def compute_and_align_subject_gradients(matrices_dict, group_gradient, labels, lh_idx, rh_idx, n_components=3):
    """
    For each subject in matrices_dict, compute the subject gradients and align them
    to the provided group_gradient using Procrustes alignment.
    
    Returns:
      - aligned: Dictionary mapping subject IDs to their aligned gradient matrices.
    """
    aligned = {}
    for subj, matrix in matrices_dict.items():
        subj_grad = compute_subject_gradients(matrix, labels, lh_idx, rh_idx, n_components)
        # Align using BrainSpace's procrustes_alignment
        aligned_grad = procrustes_alignment([subj_grad], group_gradient)[0]

        # Scale each gradient (i.e., each column) between 0 and 1
        for i in range(aligned_grad.shape[1]):
            col = aligned_grad[:, i]
            # Avoid division by zero if the column has constant values
            if np.ptp(col) != 0:
                aligned_grad[:, i] = (col - np.min(col)) / np.ptp(col)
            else:
                aligned_grad[:, i] = 0.0  # or leave as is
        aligned[subj] = aligned_grad
        print(f"Aligned gradients computed for {subj}")
    
    return aligned

# Compute aligned gradients for patients and older healthy controls
aligned_patients_gradients = compute_and_align_subject_gradients(
    patients_matrices, group_gradients_combined, parcel_labels, lh_indices, rh_indices, n_components=3
)
aligned_controls_gradients = compute_and_align_subject_gradients(
    older_hc_matrices, group_gradients_combined, parcel_labels, lh_indices, rh_indices, n_components=3
)


















# Spin permutation on gradients [one properly computed...]

from brainspace.null_models import SpinPermutations
from brainspace.datasets import load_conte69
from scipy.stats import spearmanr
# Load sphere surfaces for spin testing
sphere_lh, sphere_rh = load_conte69(as_sphere=True)

# Ensure the shapes of the sphere surfaces match the gradient mappings
n_vertices_lh = sphere_lh.points.shape[0]
n_vertices_rh = sphere_rh.points.shape[0]

# Map the averaged gradients to surfaces
patient_gradients_mapped = map_to_labels(patient_gradients_combined[:, 2], labeling, mask=labeling != 0)
older_hc_gradients_mapped = map_to_labels(older_hc_gradients_combined[:, 2], labeling, mask=labeling != 0)

# Ensure gradient mapping matches surface resolution
assert len(patient_gradients_mapped) == n_vertices_lh + n_vertices_rh, "Gradient mapping size mismatch."

# Spin Permutations
n_rand = 5000  # Number of random permutations
sp = SpinPermutations(n_rep=n_rand, random_state=0)
sp.fit(sphere_lh, points_rh=sphere_rh)  # Fit to left and right hemisphere spheres

# Rotate the PD gradient mappings
patients_rotated = np.hstack(
    sp.randomize(patient_gradients_mapped[:n_vertices_lh], patient_gradients_mapped[n_vertices_lh:])
)

# Compute observed correlation between gradients
mask = ~np.isnan(older_hc_gradients_mapped) & ~np.isnan(patient_gradients_mapped)
r_obs, pv_obs = spearmanr(older_hc_gradients_mapped[mask], patient_gradients_mapped[mask])

# Compute correlations for permuted gradients
r_spin = np.empty(n_rand)
for i, perm in enumerate(patients_rotated):
    mask_rot = mask & ~np.isnan(perm)
    r_spin[i] = spearmanr(older_hc_gradients_mapped[mask_rot], perm[mask_rot])[0]

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
    array_name=[patient_gradients_mapped, older_hc_gradients_mapped],
    size=(1200, 600),
    cmap='RdBu',
    color_bar=True,
    label_text=['Patients Tertiary Gradient', 'Older HC Tertiary Gradient'],
    zoom=1.05
)
plt.show()



# I also want to calculate regional and network measures of eccentricity.









