#! /bin/bash

export FREESURFER_HOME=/Applications/freesurfer
source $FREESURFER_HOME/SetUpFreeSurfer.sh

# Compute FOD response functions for DWI
base=~/imaging/datasets/mrgfus
raw=${base}/rawdata
dwi=${base}/derivatives/dwi
fba=${base}/derivatives/fba
sf=${base}/derivatives/study_files
# bash ~/Desktop/flip.sh
hemis=${sf}/fba/hemis.txt # sub-id hemi

for dir in ${fba}/data/flip/sub-*; do
  sub=$(basename ${dir})

  for ses in ses-01 ses-02 ses-03; do
    eddy=${dwi}/${sub}/${ses}/eddy
    
    fod=${fba}/data/flip/${sub}/${ses}/fod
    
    if [ ! -d ${eddy} ]; then
      echo "Skipping subject ${sub} ${ses} - this directory is not found!"
      continue
    fi
  
    echo "Processing subject ${sub} for session ${ses}"
    mkdir -p ${fod}

    # Flip x-axis
    mrtransform \
      ${eddy}/${sub}_dwi_edc.nii.gz \
      ${fod}/tmp_${sub}_dwi.mif \
      -flip 0 \
      -fslgrad \
      ${eddy}/${sub}_dwi_edc.eddy_rotated_bvecs \
      ${dwi}/${sub}/${ses}/${sub}_${ses}_acq-dwi.bval
  
    # Bias correct
    dwibiascorrect \
      ants \
      ${fod}/tmp_${sub}_dwi.mif \
      ${fod}/tmp_${sub}_dwi_dn.mif

    # Compute average tissue-response function
    dwi2response \
      dhollander \
      ${fod}/tmp_${sub}_dwi_dn.mif \
      ${fod}/response_wm.txt \
      ${fod}/response_gm.txt \
      ${fod}/response_csf.txt

    # Regrid to voxel size 1.25
    mrgrid \
      ${fod}/tmp_${sub}_dwi_dn.mif \
      regrid -vox 1.25 \
      ${fod}/${sub}_dwi_us.mif

    # Extract upsampled b0
    dwiextract \
      ${fod}/${sub}_dwi_us.mif \
      - -bzero | mrmath - mean ${fod}/tmp_${sub}_b0_us.mif \
      -axis 3

    mrconvert \
      ${fod}/tmp_${sub}_b0_us.mif \
      ${fod}/tmp_${sub}_b0_us.nii.gz

    # Brain extract and create mask
    mri_synthstrip \
      -i ${fod}/tmp_${sub}_b0_us.nii.gz \
      -o ${fod}/tmp_${sub}_b0_brain_us.nii.gz \
      -m ${fod}/tmp_${sub}_b0_brain_mask_us.nii.gz \
      -t 12

    # Convert brain mask to mif format
    mrconvert \
      ${fod}/tmp_${sub}_b0_brain_mask_us.nii.gz \
      ${fod}/${sub}_b0_brain_mask_us.mif \
      -force

    # Create individual tissue response functions
    source /Users/a9ws/mrtrix_env/bin/activate
    ss3t_csd_beta1 \
      ${fod}/${sub}_dwi_us.mif \
      ${fba}/rfuncs/group_response_wm.txt \
      ${fod}/tmp_${sub}_wmfod.mif \
      ${fba}/rfuncs/group_response_gm.txt \
      ${fod}/tmp_${sub}_gmfod.mif \
      ${fba}/rfuncs/group_response_csf.txt \
      ${fod}/tmp_${sub}_csffod.mif \
      -mask ${fod}/${sub}_b0_brain_mask_us.mif
    deactivate

    mtnormalise \
      ${fod}/tmp_${sub}_wmfod.mif ${fod}/${sub}_wmfod.mif \
      ${fod}/tmp_${sub}_gmfod.mif ${fod}/${sub}_gmfod.mif \
      ${fod}/tmp_${sub}_csffod.mif ${fod}/${sub}_csffod.mif \
      -mask ${fod}/${sub}_b0_brain_mask_us.mif

    # Clean up
    rm \
      ${fod}/${sub}_dwi_us.mif \
      rm ${fod}/tmp*

    itemp=${fba}/template/intra-temps/${sub}
    rm -r ${itemp}
    mkdir -p ${itemp}/fods ${itemp}/masks

    fod=${fba}/data/flip/${sub}/${ses}/fod
    
    fod_in=${fod}/${sub}_wmfod.mif
    fod_out=${itemp}/fods/${sub}_${ses}.mif
    mask_in=${fod}/${sub}_b0_brain_mask_us.mif
    mask_out=${itemp}/masks/${sub}_${ses}.mif

    if [ -f ${fod_in} ]; then
      echo "  Found FOD for ${ses}, linking to intra-temp folder."
      ln -sf ${fod_in} ${fod_out}
      ln -sf ${mask_in} ${mask_out}
    else
      echo "  No FOD found for ${ses}, skipping."
    fi
  done

  # Run population_template with rigid registration
  # Was -type rigid , now running rigid_affine_nonlinear, as default
  population_template \
    ${itemp}/fods/ \
    ${itemp}/fods/${sub}_itemp.mif \
    -mask_dir ${itemp}/masks/ \
    -voxel_size 1.25 \
    -warp_dir ${itemp}/xfms \
    -force
done


# WILL THEN NEED TO RUN THESE THROUGH INTRA_TEMPS ETC

