#! /bin/bash

der=/Users/neuero-239/imaging/datasets/mrgfus/derivatives
# Subcortical volumes

for ses in ses-01 ses-02 ses-03; do
  # Use an array to store matching directories
  subjects=(${der}/data/anat/sub-*)
  for dir in "${subjects[@]}"; do
      sub=$(basename "$dir")
      mkdir -p ${der}/projects/FIRST/${sub}/${ses}
      
      run_first_all -d -b -i ${dir}/${ses}/${sub}_${ses}_acq-T1w_brain.nii.gz \
      -o ${der}/projects/FIRST/${sub}/${ses}
  done
done
