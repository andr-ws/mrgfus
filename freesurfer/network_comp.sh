#!/bin/bash

# Code description:
# Computes .annot files for Schaefer parcellation for each subject.
# Computes Freesurfer to native xfm to get the parcellation in native space.

raw=~/imaging/datasets/mrgfus/rawdata
der=~/imaging/datasets/mrgfus/derivatives
dwi=${der}/dwi
anat=${der}/anat

export SUBJECTS_DIR=${der}/freesurfer

# LUT path
lut_path=< /point/2/here >

for dir in ${dwi}/sub-*; do
  sub=$(basename ${dir})

  for ses in ses-preop ses1-postop ses2-postop; do
    ses_dir=${SUBJECTS_DIR}/${sub}/${ses}
    mkdir ${dir}/${ses}/fba/tractogram
    
    # Resample the Schaefer parcellation to subject surface - generate annot files for multiple resolutions
    for rois in 400 600 800; do
      for h in lh rh; do
        mri_surf2surf \
          --hemi ${h} \
          --srcsubject fsaverage \
          --trgsubject ${sub}_${ses} \
          --sval-annot ${SUBJECTS_DIR}/fsaverage/label/${h}.Schaefer2018_${rois}Parcels_7Networks_order.annot \
          --tval ${ses_dir}/label/${h}.Schaefer2018_${rois}Parcels_7Networks_order.annot
      done
      
      mri_aparc2aseg \
        --s ${sub}_${ses} \
        --annot Schaefer2018_${rois}Parcels_7Networks_order \
        ${ses_dir}/mri/Schaefer2018_${rois}Parcels_7Networks_order.mgz
      
      mrconvert \
        ${ses_dir}/mri/Schaefer2018_${rois}Parcels_7Networks_order.mgz \
        ${dir}/${ses}/fba/tractogram/Schaefer2018_${rois}Parcels_7Networks_order.mif

      labelconvert \
        ${dir}/${ses}/fba/tractogram/Schaefer2018_${rois}Parcels_7Networks_order.mif \
        ${lut_path} \
        ${dir}/${ses}/fba/tractogram/Schaefer2018_${rois}Parcels_7Networks_order.txt \
        ${dir}/${ses}/fba/tractogram/Schaefer2018_${rois}Parcels_7Networks_order_conv.mif

      # MIND network output
      mkdir ${SUBJECTS_DIR}/${sub}/${ses}/MIND
        
      # Compute MIND networks from freesurfer output
      # Assumes MIND code has been sourced [git clone https://github.com/isebenius/MIND.git ]
      python3 ~/imaging/code/projects/mrgfus/MIND_exec.py ${SUBJECTS_DIR}/${sub} ${ses} ${rois}
    
    done

    # Use brain extracted FS output for ACT
    5ttgen \
      fsl \
      ${SUBJECTS_DIR}/${ses}/${sub}_brain.nii.gz \
      ${dir}/${ses}/fba/tractogram/${sub}_5tt.mif \
      -premasked \
      -force

    5tt2gmwmi \
      ${dir}/${ses}/fba/tractogram/${sub}_5tt.mif \
      ${dir}/${ses}/fba/tractogram/${sub}_gmwmi.mif \
      -force

    tckgen \
      -angle 22.5 \
      -maxlen 250 \
      -minlen 10 \
      -power 1.0 \
      ${dir}/${ses}/fba/fod/${sub}_wmfod.mif \
      -seed_gmwmi ${dir}/${ses}/fba/tractogram/${sub}_gmwmi.mif \
      -cutoff 0.10  \
      -select 20000000 \
      ${dir}/${ses}/fba/tractogram/${sub}.tck \
      -act ${dir}/${ses}/fba/tractogram/${sub}_5tt.mif \
      -force

    tcksift2 \
      ${dir}/${ses}/fba/tractogram/${sub}_SIFT.tck \
      ${dir}/${ses}/fba/fod/${sub}_wmfod.mif \
      ${dir}/${ses}/fba/tractogram/${sub}_tck-weights.txt \
      -act ${dir}/${ses}/fba/tractogram/${sub}_5tt.mif \
      -remove_untracked

    tck2connectome \
      ${dir}/${ses}/fba/tractogram/${sub}_SIFT.tck \
      parc_image.mif \
      ${dir}/${ses}/fba/tractogram/${sub}_Schaefer2018_400Parcels_7Networks.csv \
      -symmetric \
      -tck_weights_in ${dir}/${ses}/fba/tractogram/${sub}_tck-weights.txt
    # Remove original tractogram
    #rm ${dir}/fba/tractograms/${sub}.tck
  done
done
