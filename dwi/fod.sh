#! /bin/bash

export FREESURFER_HOME=/Applications/freesurfer
source $FREESURFER_HOME/SetUpFreeSurfer.sh

# Compute FOD response functions for DWI

base=~/imaging/datasets/mrgfus
raw=${base}/rawdata
dwi=${base}/derivatives/dwi

for dir in ${dwi}/sub-*; do
  sub=$(basename ${dir})

  if [ ! -d "${eddy}" ]; then
    echo "Skipping subject ${sub}, dwi/eddy directory not found."
    continue
  fi

  for ses in ses-preop ses1-postop ses2-postop; do
  
    echo "Processing subject ${sub} for session ${ses}"
    bval=${raw}/${sub}/${ses}/dwi/${sub}_dwi.bval
    mkdir -p ${dir}/${ses}/fod
  
    # Convert to mif format
    mrconvert \
      ${dir}/${ses}/eddy/${sub}_dwi_edc.nii.gz \
      ${dir}/${ses}/fod/${sub}_dwi.mif \
      -fslgrad \
      ${eddy}/${sub}_dwi_edc.eddy_rotated_bvecs \
      ${bval} \
      -force
    
    # Bias field correct
    dwibiascorrect \
      ants \
      ${dir}/${ses}/fod/${sub}_dwi.mif \
      ${dir}/${ses}/fod/tmp_${sub}_dwi_dn.mif -force

    # Compute average tissue-response function
    dwi2response \
      dhollander \
      ${dir}/${ses}/fod/tmp_${sub}_dwi_dn.mif \
      ${dir}/${ses}/fod/response_wm.txt \
      ${dir}/${ses}/fod/response_gm.txt \
      ${dir}/${ses}/fod/response_csf.txt -force

    # Regrid to voxel size 1.25
    mrgrid \
      ${dir}/${ses}/fod/tmp_${sub}_dwi_dn.mif \
      regrid -vox 1.25 \
      ${dir}/${ses}/fod/tmp_${sub}_dwi_us.mif -force

    # Extract upsampled b0 for ACT
    dwiextract \
      ${dir}/${ses}/fod/tmp_${sub}_dwi_us.mif \
      - -bzero | mrmath - mean ${dir}/${ses}/fod/${sub}_b0_us.mif \
      -axis 3

    mrconvert ${dir}/${ses}/fod/${sub}_b0_us.mif ${dir}/${ses}/fod/${sub}_b0_us.nii.gz

    # Brain extract and create mask
    mri_synthstrip \
      -i ${dir}/${ses}/fod/${sub}_b0_us.nii.gz \
      -o ${dir}/${ses}/fod/${sub}_b0_brain_us.nii.gz \
      -m ${dir}/${ses}/fod/${sub}_b0_brain_mask_us.nii.gz

    # Convert brain mask to mif format
    mrconvert \
      ${dir}/${ses}/fod/${sub}_b0_brain_mask_us.nii.gz \
      ${dir}/${ses}/fod/${sub}_b0_brain_mask_us.mif

    rm ${dir}/${ses}/fod/*tmp*
      
  done
done

fba=${base}/derivatives/fba
mkdir -p ${fba}/response

# Create tissue RF's

for tissue in wm gm csf; do
  responsemean ${dwi}/*/*/fod/response_${tissue}.txt ${fba}/response/group_response_${tissue}.txt
done

for dir in ${dwi}/sub-*; do
  sub=$(basename ${dir})

  for ses in ses-preop ses1-postop ses2-postop; do
    ss3t_csd_beta1 \
      ${dir}/${ses}/fod/tmp_${sub}_dwi_us.mif \
      ${fba}/response/group_response_wm.txt \
      ${dir}/${ses}/fod/tmp_${sub}_wmfod.mif \
      ${fba}/response/group_response_gm.txt \
      ${dir}/${ses}/fod/tmp_${sub}_gmfod.mif \
      ${fba}/response/group_response_csf.txt \
      ${dir}/${ses}/fod/tmp_${sub}_csffod.txt \
      -mask ${dir}/${ses}/fod/${sub}_b0_brain_mask_us.mif

    mtnormalise \
      ${dir}/${ses}/fod/tmp_${sub}_wmfod.mif ${dir}/${ses}/fod/${sub}_wmfod.mif \
      ${dir}/${ses}/fod/tmp_${sub}_gmfod.mif ${dir}/${ses}/fod/${sub}_gmfod.mif
      ${dir}/${ses}/fod/tmp_${sub}_csffod.mif ${dir}/${ses}/fod/${sub}_csffod.mif

    rm ${dir}/${ses}/fod/*tmp*
  done
done

