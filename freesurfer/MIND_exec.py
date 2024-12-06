import sys
import os
from MIND import compute_MIND

if len(sys.argv) < 4:
    print("Usage: python MIND_exec.py <path_to_surf_dir> <session> <rois>")
    sys.exit(1)

# Get arguments
path_to_surf_dir = sys.argv[1]
session = sys.argv[2]
rois = sys.argv[3]

# Validate path_to_surf_dir
if not os.path.exists(path_to_surf_dir):
    print(f"Error: Directory {path_to_surf_dir} does not exist.")
    sys.exit(1)

# Features and parcellation
features = ['CT', 'MC', 'Vol', 'SD', 'SA']
parcellation = f"Schaefer2018_{rois}Parcels_7Networks_order.annot"

# Compute MIND network
print(f"Computing MIND network for session: {session}, resolution: {rois}")
MIND = compute_MIND(path_to_surf_dir, features, parcellation)

# Prepare output directory
output_dir = os.path.join(path_to_surf_dir, session, "MIND")
os.makedirs(output_dir, exist_ok=True)

# Save MIND network to file
output_file = os.path.join(output_dir, f"MIND_Schaefer2018_{rois}Parcels_7Networks_{session}.csv")
MIND.to_csv(output_file)
print(f"MIND network saved to {output_file}")
