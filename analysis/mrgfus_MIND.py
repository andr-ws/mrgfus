import os
import sys
import pandas as pd
from MIND import compute_MIND

# Set the path to the FreeSurfer subjects directory
freesurfer_dir = "/Users/neuero-239/imaging/datasets/mrgfus/derivatives/data/freesurfer"

# Define MIND features and parcellation
features = ['CT', 'MC', 'Vol', 'SD', 'SA']
parcellation = 'Schaefer2018_400Parcels_7Networks_order'

# Define sessions to iterate over
sessions = ['ses-01', 'ses-02', 'ses-03']

# Output directory for results
output_dir = "/Users/neuero-239/imaging/datasets/mrgfus/derivatives/mind_networks"
os.makedirs(output_dir, exist_ok=True)

# Get list of subjects
subjects = [sub for sub in os.listdir(freesurfer_dir) if os.path.isdir(os.path.join(freesurfer_dir, sub))]

# Loop over subjects and sessions
for subject in subjects:
    for session in sessions:
        # Define session path
        session_path = os.path.join(freesurfer_dir, subject, session)
        
        # Check if session directory exists before processing
        if not os.path.exists(session_path):
            print(f"Skipping {subject} {session}: Directory not found.")
            continue
        
        try:
            print(f"Processing {subject}, {session}...")

            # Compute MIND network
            mind_network = compute_MIND(session_path, features, parcellation)

            # Save the output as CSV
            output_file = os.path.join(output_dir, f"{subject}_{session}_MIND.csv")
            mind_network.to_csv(output_file)
            print(f"Saved MIND network for {subject}, {session} at {output_file}")

        except Exception as e:
            print(f"Error processing {subject}, {session}: {e}")

print("Batch processing complete!")
