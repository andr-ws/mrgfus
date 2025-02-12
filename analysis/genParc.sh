#!/bin/bash

# Code description:
# Computes .annot files for Schaefer parcellations for each subject and puts them in upsampled dwi space

der=~/imaging/datasets/mrgfus/derivatives
export SUBJECTS_DIR=${der}/data/freesurfer

for dir in ${der}/data/freesurfer/sub-*; do
  sub=$(basename ${dir})
  parc=${dir}/parcellations
  
  mkdir ${parc}

  # Resample the Schaefer parcellations to subject surface - generate annot files for multiple resolutions
  # UPDATE: resample parcellations to synthetic b0 surface
  for rois in 400; do
    for h in lh rh; do
        mri_surf2surf \
        --hemi ${h} \
        --srcsubject fsaverage \
        --trgsubject ${dwi}/freesurfer \
        --sval-annot ${SUBJECTS_DIR}/fsaverage/label/${h}.Schaefer2018_${rois}Parcels_7Networks_order.annot \
        --tval ${dwi}/freesurfer/label/${h}.Schaefer2018_${rois}Parcels_7Networks_order.annot
    done
          
    # Convert annots to volume
    mri_aparc2aseg --s ${dwi}/freesurfer --o ${dwi}/parcellations/${sub}_Schaefer2018_${rois}Parcels_7Networks.nii.gz --annot Schaefer2018_${rois}Parcels_7Networks_order

    # Reorient the parcellations
    fslreorient2std ${dwi}/parcellations/${sub}_Schaefer2018_${rois}Parcels_7Networks.nii.gz ${dwi}/parcellations/${sub}_Schaefer2018_${rois}Parcels_7Networks.nii.gz

    # bval and bvec files
    bval=${dwi}/${sub}_dwi.bval
    bvec=${dwi}/eddy/${sub}_dwi_edc.eddy_rotated_bvecs

    # Regrid dwi to 1.25 iso
    mrgrid ${dwi}/eddy/${sub}_dwi_edc.nii.gz regrid -vox 1.25 ${dwi}/parcellations/tmp_data.nii.gz -force
      
    # Extract an upsampled b0
    dwiextract ${dwi}/parcellations/tmp_data.nii.gz - -bzero -fslgrad ${bvec} ${bval} | mrmath - mean ${dwi}/parcellations/tmp_b0_us.nii.gz -axis 3 -force
     
    # Brain extract and mask the upsampled b0
    bet2 ${dwi}/parcellations/tmp_b0_us.nii.gz ${dwi}/parcellations/b0_us_brain -m
    
    rm ${dwi}/parcellations/tmp*
    
    # Compute (rigid) us_b0 to synthetic space xfm
    mri_convert ${dwi}/freesurfer/mri/synthSR.norm.mgz ${dwi}/parcellations/sb0_freesurfer.nii.gz
    fslreorient2std ${dwi}/parcellations/sb0_freesurfer.nii.gz ${dwi}/parcellations/sb0_freesurfer.nii.gz

    antsRegistrationSyN.sh -d 3 -t r -n 16 \
    -f ${dwi}/parcellations/sb0_freesurfer.nii.gz \
    -m ${dwi}/parcellations/b0_us_brain.nii.gz \
    -o ${dwi}/parcellations/b0_to_sb0_

    # Transform parcellations to dwi space
    antsApplyTransforms -d 3 \
    -r ${dwi}/parcellations/b0_us_brain.nii.gz \
    -i ${dwi}/parcellations/${sub}_Schaefer2018_${rois}Parcels_7Networks.nii.gz \
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
