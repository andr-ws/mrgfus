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
    mkdir -p ${dir}/${ses}/fixels
  
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
fod2fixel \
-mask ${fba}/template/template_mask.mif \
-fmls_peak_value 0.1 \                         #<- correct value?
${fba}/template/wmfod_template.mif \
${fba}/template/fixel_mask

# Segment FOD images to estimate fixels and their AFD
for dir in ${fba}/data/sub-*; do
  sub=$(basename ${dir})
  for ses in ses-01 ses-02 ses-03; do
  
    fod2fixel \
    -mask ${fba}/template/template_mask.mif \
    ${dir}/${ses}/fod/${sub}_fod-template_NOT_REORIENTED.mif \
    ${dir}/${ses}/fixels/${sub}_fixel-template_NOT_REORIENTED \
    -afd ${ses}_fd.mif \            #<- not sure where this goes?
    -force
   
    # Reorient fixels
    fixelreorient \
    ${dir}/${ses}/fixels/${sub}_fixel-template_NOT_REORIENTED \
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
