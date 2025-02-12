#!/bin/bash

# Code description:
# Computes .annot files for Schaefer parcellations for each subject and puts them in upsampled dwi space

der=~/imaging/datasets/mrgfus/derivatives
export SUBJECTS_DIR=${der}/data/freesurfer

for dir in ${der}/data/freesurfer/sub-*; do
  sub=$(basename ${dir})
  for ses in ses-01 ses-02 ses-03
    parc=${dir}/${ses}/parcellations
  
    mkdir ${parc}

    # Resample the Schaefer parcellations to subject surface - generate annot files for multiple resolutions (done in SchaeferResample.sh)      
    # Convert annots to volume
    mri_aparc2aseg --s ${dir}/${ses} --o ${parc}/${sub}_Schaefer2018_400Parcels_7Networks.nii.gz --annot Schaefer2018_400Parcels_7Networks_order

    # Reorient the parcellations
    fslreorient2std ${parc}/${sub}_Schaefer2018_400Parcels_7Networks.nii.gz ${parc}/${sub}_Schaefer2018_400Parcels_7Networks.nii.gz

    # bval and bvec files
    bval=${der}/data/dwi/${sub}/${ses}/${sub}_${ses}_acq-dwi.bval
    bvec=${der}/data/dwi/${sub}/${ses}/eddy/${sub}_dwi_edc.eddy_rotated_bvecs
      
    # Extract an upsampled b0
    dwiextract ${der}/projects/fba/data/${sub}/${ses}/fod/${sub}_dwi_us.mif - -bzero -fslgrad ${bvec} ${bval} | mrmath - mean ${parc}/tmp_${sub}_b0_us.nii.gz -axis 3 -force
     
    # Brain extract and mask the upsampled b0
    bet2 ${parc}/tmp_${sub}_b0_us.nii.gz ${parc}/${sub}_b0_us_brain -m
    
    # Compute (rigid) us_b0 to synthetic space xfm
    mri_convert ${dir}/${ses}/mri/brain.norm.mgz ${parc}/tmp_${sub}_fs_brain.nii.gz
    fslreorient2std ${parc}/tmp_${sub}_fs_brain.nii.gz ${parc}/${sub}_fs_brain.nii.gz

    rm ${parc}/tmp*

    antsRegistrationSyN.sh -d 3 -t r -n 16 \
    -f ${parc}/${sub}_fs_brain.nii.gz \
    -m ${parc}/${sub}_b0_us_brain.nii.gz \
    -o ${parc}/${sub}_b0_to_fs_

    # Transform parcellations to dwi space
    antsApplyTransforms -d 3 \
    -r ${parc}/${sub}_b0_us_brain.nii.gz \
    -i ${parc}/${sub}_Schaefer2018_${rois}Parcels_7Networks.nii.gz \
    -o ${dwi}/parcellations/${sub}_Schaefer2018_${rois}Parcels_7Networks_space-dwi.nii.gz \
    -t [${dwi}/parcellations/b0_to_sb0_0GenericAffine.mat,1] \
    -n GenericLabel[linear]

    # Remove non-cortical labels
    fslmaths \
    ${dwi}/parcellations/${sub}_Schaefer2018_${rois}Parcels_7Networks_space-dwi.nii.gz \
    -thr 1000 \
    ${dwi}/parcellations/${sub}_Schaefer2018_${rois}Parcels_7Networks_space-dwi.nii.gz

    # Remove medial walls
    fslmaths ${dwi}/parcellations/${sub}_Schaefer2018_${rois}Parcels_7Networks_space-dwi.nii.gz \
    -thr 1000 -uthr 1000 ${dwi}/parcellations/tmp_Schaefer2018_${rois}Parcels_7Networks_space-dwi_exclude1.nii.gz
    fslmaths ${dwi}/parcellations/${sub}_Schaefer2018_${rois}Parcels_7Networks_space-dwi.nii.gz \
    -thr 2000 -uthr 2000 ${dwi}/parcellations/tmp_Schaefer2018_${rois}Parcels_7Networks_space-dwi_exclude2.nii.gz
    fslmaths ${dwi}/parcellations/tmp_Schaefer2018_${rois}Parcels_7Networks_space-dwi_exclude1.nii.gz \
    -add ${dwi}/parcellations/tmp_Schaefer2018_${rois}Parcels_7Networks_space-dwi_exclude2.nii.gz \
    ${dwi}/parcellations/tmp_Schaefer2018_${rois}Parcels_7Networks_space-dwi_exclude.nii.gz
    fslmaths ${dwi}/parcellations/${sub}_Schaefer2018_${rois}Parcels_7Networks_space-dwi.nii.gz \
    -sub ${dwi}/parcellations/tmp_Schaefer2018_${rois}Parcels_7Networks_space-dwi_exclude.nii.gz \
    ${dwi}/parcellations/${sub}_Schaefer2018_${rois}Parcels_7Networks_space-dwi.nii.gz

    # Clean up
    rm ${dwi}/parcellations/tmp*.nii.gz
    
  done
done
