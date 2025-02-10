#!/bin/bash

# Code description:
# Computes .annot files for Schaefer parcellations for each subject

der=~/imaging/datasets/mrgfus/derivatives
export SUBJECTS_DIR=${der}/data/freesurfer
ln -s /Applications/freesurfer/dev/subjects/fsaverage /Users/neuero-239/imaging/datasets/mrgfus/derivatives/data/freesurfer/

for dir in ${SUBJECTS_DIR}/sub-*; do
  sub=$(basename ${dir})
  for h in lh rh; do
    mri_surf2surf \
    --hemi ${h} \
    --srcsubject fsaverage \
    --trgsubject ${dir} \
    --sval-annot ${SUBJECTS_DIR}/fsaverage/label/${h}.Schaefer2018_${rois}Parcels_7Networks_order.annot \
    --tval ${dir}/label/${h}.Schaefer2018_${rois}Parcels_7Networks_order.annot
  done
done
