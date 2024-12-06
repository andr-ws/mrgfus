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

  for ses in ses-preop ses1-postop ses2-postop; do
    # Resample the Schaefer parcellation to subject surface - generate annot files for multiple resolutions
    for rois in 400 600 800; do
      for h in lh rh; do
        mri_surf2surf \
          --hemi ${h} \
          --srcsubject fsaverage \
          --trgsubject ${sub}_${ses} \
          --sval-annot ${SUBJECTS_DIR}/fsaverage/label/${h}.Schaefer2018_${rois}Parcels_7Networks_order.annot \
          --tval ${SUBJECTS_DIR}/${sub}/label/${h}.Schaefer2018_${rois}Parcels_7Networks_order.annot
      done
      
          
    # Extrapolate FS metrics (MIND networks)
    --> Execute code here.

    # ACT
    mkdir ${dir}/fba/tractograms

    # NOT SURE IF REQUIRED YET
    # Map b0 to T1 (rigid)
    #antsRegistrationSyN.sh \
      #  -d 3 \
      #  -f T1_bias-cor \
      #  -m eddy_b0 \
      #  -t r
    # Invert
    #antsApplyTransforms \
      #  -i T1_bias-cor \
      #  -r eddy_b0 \
    
    # Use brain extracted FS output
    5ttgen \
      fsl \
      ${SUBJECTS_DIR}/${ses}/${sub}_brain.nii.gz \
      ${dir}/fba/tractograms/${sub}_5tt.mif \
      -premasked \
      -force

    5tt2gmwmi \
      ${dir}/fba/tractograms/${sub}_5tt.mif \
      ${dir}/fba/tractograms/${sub}_gmwmi.mif \
      -force

    tckgen \
      -angle 22.5 \
      -maxlen 250 \
      -minlen 10 \
      -power 1.0 \
      ${dir}/fba/proc/${sub}_wmfod_norm.mif \
      -seed_gmwmi ${dir}/fba/proc/${sub}_gmwmi.mif \
      -cutoff 0.10  \
      -select 20000000 \
      ${dir}/fba/tractograms/${sub}.tck \
      -act ${dir}/tractograms/${sub}_5tt.mif \
      -force

  # SIFT tractogram
  #tcksift \
  #  ${dir}/fba/tractograms/${sub}.tck \
  #  ${dir}/fba/proc/${sub}_wmfod_norm.mif \
  #  ${dir}/fba/tractograms/${sub}_SIFT.tck \
  #  -term_number 2000000 \
  #  -force

    tcksift2 \
      ${dir}/fba/tractograms/${sub}_SIFT.tck \
      ${dir}/fba/proc/${sub}_wmfod_norm.mif \
      ${dir}/fba/tractograms/${sub}.txt \
      -act ${dir}/tractograms/${sub}_5tt.mif

    tck2connectome \
      ${dir}/fba/tractograms//ACT/${sub}_SIFT.tck \
      parc_image.mif \
      ${fba}/subjects/${sub}/ACT/Schaefer_tck_400.csv \
      -symmetric \
      -tck_weights_in ${dir}/fba/tractograms/${sub}.txt

    # Remove original tractogram
    #rm ${dir}/fba/tractograms/${sub}.tck
  done
done
