#!/bin/bash

# Code description:
# Computes .annot files for Schaefer parcellations for each subject

der=~/imaging/datasets/mrgfus/derivatives
export SUBJECTS_DIR=${der}/data/freesurfer
ln -s /Applications/freesurfer/dev/subjects/fsaverage /Users/neuero-239/imaging/datasets/mrgfus/derivatives/data/freesurfer/

for dir in ${SUBJECTS_DIR}/sub-*; do
  sub=$(basename ${dir})
  for ses in ses-01 ses-02 ses-03; do
    for h in lh rh; do
      mri_surf2surf \
      --hemi ${h} \
      --srcsubject fsaverage \
      --trgsubject ${sub}/${ses} \
      --sval-annot ${SUBJECTS_DIR}/fsaverage/label/${h}.Schaefer2018_400Parcels_7Networks_order.annot \
      --tval ${dir}/${ses}/label/${h}.Schaefer2018_400Parcels_7Networks_order.annot
    done
  done
done
