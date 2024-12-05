#!/bin/bash

# Code description:
# Computes .annot files for Schaefer parcellation for each subject.
# Computes Freesurfer to native xfm to get the parcellation in native space.

raw=~/imaging/datasets/mrgfus/rawdata
der=~/imaging/datasets/mrgfus/derivatives
dwi=${der}/dwi
anat=${der}/anat
export SUBJECTS_DIR=${der}/freesurfer

for dir in ${dwi}/sub-*; do
  sub=$(basename ${dir})
        
  # Resample the Schaefer parcellation to subject surface - generate annot files for multiple resolutions
  for rois in 400 600 800; do
    for h in lh rh; do
      mri_surf2surf \
        --hemi ${h} \
        --srcsubject fsaverage \
        --trgsubject ${sub} \
        --sval-annot ${SUBJECTS_DIR}/fsaverage/label/${h}.Schaefer2018_${rois}Parcels_7Networks_order.annot \
        --tval ${SUBJECTS_DIR}/${sub}/label/${h}.Schaefer2018_${rois}Parcels_7Networks_order.annot
    done
          
    # Convert annot to volume
    mri_aparc2aseg \
      --s ${sub} \
      --o ${SUBJECTS_DIR}/${sub}/mri/Schaefer2018_${rois}Parcels_7Networks.mgz \
      --annot Schaefer2018_${rois}Parcels_7Networks_order

    # Convert and reorient the parcellations and brain file
    mri_convert \
      --in_orientation LIA \
      --out_orientation RAS \
      -it mgz \
      -ot nii \
      ${SUBJECTS_DIR}/${sub}/mri/Schaefer2018_${rois}Parcels_7Networks.mgz \
      ${SUBJECTS_DIR}/${sub}/mri/Schaefer2018_${rois}Parcels_7Networks.nii.gz

    # Apply native transform to the Schaefer parcellation
    mkdir ${anat}/${sub}/parc
      
    antsApplyTransforms \
      -d 3 \
      -i ${SUBJECTS_DIR}/${sub}/mri/Schaefer2018_${rois}Parcels_7Networks.nii.gz \
      -r ${anat}/${sub}/proc/${sub}_desc-res-1mm_bfc_T1w_brain.nii.gz \
      -o ${anat}/${sub}/parc/${sub}_Schaefer2018_${rois}Parcels_7Networks_T1w.nii.gz \
      -n NearestNeighbor \
      -t ${SUBJECTS_DIR}/${sub}/mri/transforms/fs_2_nat_0GenericAffine.mat
  done
done













