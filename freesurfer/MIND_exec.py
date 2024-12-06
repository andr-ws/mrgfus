import sys
import os
from MIND import compute_MIND

if len(sys.argv) < 3:
    print("Usage: python MIN_exec.py <path_to_surf_dir> <session>")
    sys.exit(1)

# Get arguments
path_to_surf_dir = sys.argv[1]
session = sys.argv[2]

# Features and parcellation
features = ['CT', 'MC', 'Vol', 'SD', 'SA']

# THIS NEEDS TO ITERATE THROUGH 200 400 600 PARCELS
parcellation = 'Schaefer2018_400Parcels_7Networks_order.annot'

# Compute MIND network
print(f"Computing MIND network for session: {session}")
MIND = compute_MIND(path_to_surf_dir, features, parcellation)

# Prepare output directory
output_dir = os.path.join(path_to_surf_dir, session, "MIND")
os.makedirs(output_dir, exist_ok=True)

# Save MIND network to file
output_file = os.path.join(output_dir, f"MIND_Schaefer2018_400Parcels_7Networks_{session}.csv")
MIND.to_csv(output_file)
print(f"MIND network saved to {output_file}")
