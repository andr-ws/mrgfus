#! /bin/bash
set +H  # Disable history expansion
export FREESURFER_HOME=/Applications/freesurfer
source $FREESURFER_HOME/SetUpFreeSurfer.sh

# Compute FOD response functions for DWI
base=~/imaging/datasets/mrgfus
raw=${base}/rawdata
dwi=${base}/derivatives/data/dwi
fba=${base}/derivatives/projects/fba

for dir in ${dwi}/sub-*; do
  sub=$(basename ${dir})

  for ses in ses-01 ses-02 ses-03; do

    if [ ! -d ${dir}/${ses}/eddy ]; then
      echo "Skipping subject ${sub} ${ses} - this directory is not found!"
      continue
    fi
  
    echo "Processing subject ${sub} for session ${ses}"
    
    mkdir -p ${fba}/data/${sub}/${ses}/fod
  
    # Convert eddy data
    mrconvert \
      ${dir}/${ses}/eddy/${sub}_dwi_edc.nii.gz \
      ${fba}/data/${sub}/${ses}/fod/tmp_${sub}_dwi.mif \
      -fslgrad \
      ${dir}/${ses}/eddy/${sub}_dwi_edc.eddy_rotated_bvecs \
      ${dir}/${ses}/${sub}_${ses}_acq-dwi.bval \
    
    # Bias correct
    dwibiascorrect \
      ants \
      ${fba}/data/${sub}/${ses}/fod/tmp_${sub}_dwi.mif \
      ${fba}/data/${sub}/${ses}/fod/tmp_${sub}_dwi_dn.mif

    # Compute average tissue-response function
    dwi2response \
      dhollander \
      ${fba}/data/${sub}/${ses}/fod/tmp_${sub}_dwi_dn.mif \
      ${fba}/data/${sub}/${ses}/fod/response_wm.txt \
      ${fba}/data/${sub}/${ses}/fod/response_gm.txt \
      ${fba}/data/${sub}/${ses}/fod/response_csf.txt

    # Regrid to voxel size 1.25
    mrgrid \
      ${fba}/data/${sub}/${ses}/fod/tmp_${sub}_dwi_dn.mif \
      regrid -vox 1.25 \
      ${fba}/data/${sub}/${ses}/fod/${sub}_dwi_us.mif

    # Extract upsampled b0
    dwiextract \
      ${fba}/data/${sub}/${ses}/fod/${sub}_dwi_us.mif \
      - -bzero | mrmath - mean ${fba}/data/${sub}/${ses}/fod/tmp_${sub}_b0_us.mif \
      -axis 3

    mrconvert \
      ${fba}/data/${sub}/${ses}/fod/tmp_${sub}_b0_us.mif \
      ${fba}/data/${sub}/${ses}/fod/tmp_${sub}_b0_us.nii.gz

    # Brain extract and create mask
    mri_synthstrip \
      -i ${fba}/data/${sub}/${ses}/fod/tmp_${sub}_b0_us.nii.gz \
      -o ${fba}/data/${sub}/${ses}/fod/tmp_${sub}_b0_brain_us.nii.gz \
      -m ${fba}/data/${sub}/${ses}/fod/tmp_${sub}_b0_brain_mask_us.nii.gz

    # Convert brain mask to mif format
    mrconvert \
      ${fba}/data/${sub}/${ses}/fod/tmp_${sub}_b0_brain_mask_us.nii.gz \
      ${fba}/data/${sub}/${ses}/fod/${sub}_b0_brain_mask_us.mif

    rm ${fba}/data/${sub}/${ses}/fod/tmp*
      
  done
done

# Create group tissue response functions
mkdir -p ${fba}/rfuncs

for tissue in wm gm csf; do
  responsemean ${fba}/*/*/*/fod/response_${tissue}.txt ${fba}/rfuncs/group_response_${tissue}.txt
done

# Create individual tissue response functions
for dir in ${fba}/data/sub-*; do
  sub=$(basename ${dir})

  for ses in ses-01 ses-02 ses-03; do
    ss3t_csd_beta1 \
      ${dir}/${ses}/fod/${sub}_dwi_us.mif \
      ${fba}/rfuncs/group_response_wm.txt \
      ${dir}/${ses}/fod/tmp_${sub}_wmfod.mif \
      ${fba}/rfuncs/group_response_gm.txt \
      ${dir}/${ses}/fod/tmp_${sub}_gmfod.mif \
      ${fba}/rfuncs/group_response_csf.txt \
      ${dir}/${ses}/fod/tmp_${sub}_csffod.mif \
      -mask ${dir}/${ses}/fod/${sub}_b0_brain_mask_us.mif

    mtnormalise \
      ${dir}/${ses}/fod/tmp_${sub}_wmfod.mif ${dir}/${ses}/fod/${sub}_wmfod.mif \
      ${dir}/${ses}/fod/tmp_${sub}_gmfod.mif ${dir}/${ses}/fod/${sub}_gmfod.mif \
      ${dir}/${ses}/fod/tmp_${sub}_csffod.mif ${dir}/${ses}/fod/${sub}_csffod.mif \
      -mask ${dir}/${ses}/fod/${sub}_b0_brain_mask_us.mif

    rm ${dir}/${ses}/fod/${sub}_dwi_us.mif
    rm ${dir}/${ses}/fod/tmp_${sub}_wmfod.mif ${dir}/${ses}/fod/tmp_${sub}_gmfod.mif ${dir}/${ses}/fod/tmp_${sub}_csffod.mif 
  done
done

# Create population template from subset of ses-01
template=${base}/derivatives/projects/study_files/sub_template.txt
mkdir ${fba}/template ${fba}/template/fods ${fba}/template/masks

while read -r sub; do
  ln -sf ${fba}/data/${sub}/ses-01/fod/${sub}_wmfod.mif ${fba}/template/fods/
  ln -sf ${fba}/data/${sub}/ses-01/fod/${sub}_b0_brain_mask_us.mif ${fba}/template/masks/${sub}_mask.mif
done < $template

population_template \
${fba}/template/fods/ \
-mask ${fba}/template/masks/ \
${fba}/template/wmfod_template.mif \
-voxel_size 1.25 \
-initial_alignment geometric

# Register wmFODs to template
for dir in ${fba}/data/sub-*; do
  sub=$(basename ${dir})
  
  for ses in ses-01 ses-02 ses-03; do
  
    mrregister ${dir}/${ses}/fod/${sub}_wmfod.mif -mask1 ${dir}/${ses}/fod/${sub}_b0_brain_mask_us.mif \
    ${fba}/template/wmfod_template.mif \
    -nl_warp ${dir}/${ses}/fod/${sub}-template_warp.mif \
    ${dir}/${ses}/fod/template-${sub}_warp.mif

    # Transform to template
    mrtransform ${dir}/${ses}/fod/${sub}_b0_brain_mask_us.mif \
    -warp ${dir}/${ses}/fod/${sub}-template_warp.mif \
    -interp nearest -datatype bit \
    ${dir}/${ses}/fod/${sub}_b0_mask_us-template.mif
  done
done

mrmath ${fba}/data/*/*/fod/*_b0_mask_us-template.mif \
min \
${fba}/template/template_mask.mif -datatype bit

# Define group white matter fixel mask and estimate fixel metrics for each patient
# 0.08 looks the most balanced from testing (0.06 > 0.1)
fod2fixel \
-mask ${fba}/template/template_mask.mif \
-fmls_peak_value 0.08 \
${fba}/template/wmfod_template.mif \
${fba}/template/fixel_mask

# Segment FOD images to estimate fixels and their AFD
for dir in ${fba}/data/sub-*; do
  sub=$(basename ${dir})
  for ses in ses-01 ses-02 ses-03; do

    mkdir -p ${dir}/${ses}/fixels

    mrtransform ${dir}/${ses}/fod/${sub}_wmfod.mif \
    -warp ${dir}/${ses}/fod/${sub}-template_warp.mif \
    --reorient_fod no \
    ${dir}/${ses}/fod/${sub}_wmfod_noreo.mif
    
    fod2fixel \
    -mask ${fba}/template/template_mask.mif \
    ${dir}/${ses}/fod/${sub}_wmfod_noreo.mif \
    ${dir}/${ses}/fixels/${sub}_fixel-template_noreo \
    -afd ${ses}_fd.mif \
    -force
   
    # Reorient fixels
    fixelreorient \
    ${dir}/${ses}/fixels/${sub}_fixel-template_noreo \
    ${dir}/${ses}/fod/${sub}-template_warp.mif \
    ${dir}/${ses}/fixels/${sub}_fixel-template \
    -force

    # Assign subjects fixels to template fixels
    fixelcorrespondence \
    ${dir}/${ses}/fixels/${sub}_fixel-template/fd.mif \
    ${fba}/template/fixel_mask \
    ${fba}/template/fd \
    ${sub}_${ses}.mif \
    -force

    # Compute FC metric
    warp2metric \
    ${dir}/${ses}/fod/${sub}-template_warp.mif \
    -fc ${fba}/template/fixel_mask \
    ${fba}/template/fc \
    ${sub}_${ses}.mif \
    -force
  done
done

# Copy files for, and compute log-fc
mkdir ${fba}/template/log_fc
cp ${fba}/template/fc/index.mif ${fba}/template/fc/directions.mif ${fba}/template/log_fc
for dir in ${fba}/data/sub-*; do
   sub=$(basename ${dir})
   for ses in ses-01 ses-02 ses-03; do
     mrcalc \
     ${fba}/template/fc/${sub}_${ses}.mif \
     -log ${fba}/template/log_fc/${sub}_${ses}.mif \
     -force
   done
done

# Copy files for, and compute fdc
mkdir ${fba}/template/fdc
cp ${fba}/template/fc/index.mif ${fba}/template/fc/directions.mif ${fba}/template/fdc
for dir in ${fba}/data/sub-*; do
   sub=$(basename ${dir})
   for ses in ses-01 ses-02 ses-03; do
     mrcalc \
     ${fba}/template/fd/${sub}_${ses}.mif ${fba}/template/fc/${sub}_${ses}.mif -mult \
     ${fba}/template/fdc/${sub}_${ses}.mif \
     -force
   done
done

# Compute fixel-fixel connectivity matrix
fixelconnectivity \
fixel_mask/ \
${fba}/template/tractogram_4mil_SIFT.tck matrix/ \
-force

# Smooth metric data using the fixel-fixel connectivity matrix
for metric in fd log_fc fdc; do
  fixelfilter \
  ${fba}/template/${metric} \
  smooth ${fba}/template/${metric}_smooth \
  -matrix ${fba}/template/matrix/
done

# To integrate: percent change images for postop. timepoints

#!/bin/bash

fdcdir=/Volumes/LA_4TB/mrgfus/derivatives/fba/template/fdc_smooth

# Create output directories
mkdir -p ${fdcdir}/pc/6m ${fdcdir}/pc/12m

# Loop through all ses-01 (baseline) files
for pre in ${fdcdir}/sub-*_ses-01.mif; do
  # Extract subject ID (e.g., sub-001)
  sub=$(basename $pre | cut -d '_' -f 1)

  # Define paths to sessions
  fdc_6m=${fdcdir}/${sub}_ses-02.mif
  fdc_12m=${fdcdir}/${sub}_ses-03.mif

  # Output filenames
  pc_6m=${fdcdir}/pc/6m/${sub}.mif
  pc_12m=${fdcdir}/pc/12m/${sub}.mif

  # Compute 6-month percent change if ses-02 exists
  if [ -f $fdc_6m ]; then
    echo "Computing 6m percent change for ${sub}"
    mrcalc $pre $fdc_6m -subtract $pre -divide 100 -mult $pc_6m
  else
    echo "Skipping 6m for ${sub} (missing ses-02)"
  fi

  # Compute 12-month percent change if ses-03 exists
  if [ -f $fdc_12m ]; then
    echo "Computing 12m percent change for ${sub}"
    mrcalc $pre $fdc_12m -subtract $pre -divide 100 -mult $pc_12m
  else
    echo "Skipping 12m for ${sub} (missing ses-03)"
  fi

done

# CREATE A .CSV FILE TO STORE THE MEAN VALUE PER SESSION
# IF SUBJECT DOESN'T HAVE THE SESSION, ENTRY = NA

# I want these to have rows as subject ID's and 3 columns for sessions (ses-01, ses-02, ses-03)
echo "sub-id,ses-01,ses-02,ses-03" > ${fba}/analysis/lh_hypo.csv
echo "sub-id,ses-01,ses-02,ses-03" > ${fba}/analysis/rh_hypo.csv

# Loop through each subject
for dir in ${fba}/data/sub-*; do
  sub=$(basename ${dir})
  echo "Processing $sub"

  lh_line="$sub"
  rh_line="$sub"

  # Loop through sessions
  for ses in ses-01 ses-02 ses-03; do
    fdc_mif="${fba}/template/fdc_smooth/${sub}_${ses}.mif"
    fdc_nii="${fba}/data/${sub}/${ses}/fod/${sub}_fdc_map.nii.gz"

    # Check if input file exists
    if [[ -f "$fdc_mif" ]]; then
      fixel2voxel "$fdc_mif" mean "$fdc_nii"

      for hemi in lh rh; do
        hypo_mask="${fba}/template/masks/${hemi}_hypo.nii.gz"
        fdc_masked="${fba}/data/${sub}/${ses}/fod/${sub}_fdc_${hemi}_hypo.nii.gz"

        fslmaths "$fdc_nii" -mul "$hypo_mask" "$fdc_masked"
        value=$(fslstats "$out_masked" -M)

        # Append value to the correct hemisphere line
        if [ "$hemi" == "lh" ]; then
          lh_line="$lh_line,$value"
        else
          rh_line="$rh_line,$value"
        fi
      done
    else
      # If file doesn't exist, append NA
      lh_line="$lh_line,NA"
      rh_line="$rh_line,NA"
    fi
  done

  # Append full row to CSV
  echo "$lh_line" >> ${fba}/analysis/lh_hypo.csv
  echo "$rh_line" >> ${fba}/analysis/rh_hypo.csv
done

# Hypointensity modelling

fslmaths /Volumes/LA_4TB/datasets/mrgfus/derivatives/fba/data/sub-001/ses-01/fod/sub-001_fdc_lh_hypo.nii.gz -sub /Volumes/LA_4TB/datasets/mrgfus/derivatives/fba/data/sub-001/ses-03/fod/sub-001_fdc_lh_hypo.nii.gz -div /Volumes/LA_4TB/datasets/mrgfus/derivatives/fba/data/sub-001/ses-01/fod/sub-001_fdc_lh_hypo.nii.gz -mul 100 /Volumes/LA_4TB/datasets/mrgfus/derivatives/fba/data/sub-001/ses-03/diff2.nii.gz


