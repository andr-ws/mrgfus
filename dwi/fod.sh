#! /bin/bash

export FREESURFER_HOME=/Applications/freesurfer
source $FREESURFER_HOME/SetUpFreeSurfer.sh

# Compute FOD response functions for DWI
base=~/imaging/datasets/mrgfus
raw=${base}/rawdata
dwi=${base}/derivatives/data/dwi
fba=${base}/derivatives/fba

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

# Longitudinal FBA (oh god)!

# Desc: generate rigid intra-pop templates and affine xfms for each session
# Generate a study template and write out nlin warps for those individuals
# Generate nonlinear study template warps for those not in the template
# Compose the linear and nonlinear transforms

for dir in ${fba}/data/sub-*; do
  sub=$(basename ${dir})

  mkdir -p ${fba}/template/intra-temps/${sub}/fods
  mkdir -p ${fba}/template/intra-temps/${sub}/masks

  # Loop through possible sessions
  for ses in ses-01 ses-02 ses-03; do
    fod_in=${fba}/data/${sub}/${ses}/fod/${sub}_wmfod.mif
    fod_out=${fba}/template/intra-temps/${sub}/fods/${sub}_${ses}.mif
    mask_in=${fba}/data/${sub}/${ses}/fod/${sub}_b0_brain_mask_us.mif
    mask_out=${fba}/template/intra-temps/${sub}/masks/${sub}_${ses}.mif

    if [ -f ${fod_in} ]; then
      echo "  Found FOD for ${ses}, linking to intra-temp folder."
      ln -sf ${fod_in} ${fod_out}
      ln -sf ${mask_in} ${mask_out}
    else
      echo "  No FOD found for ${ses}, skipping."
    fi
  done

  # Run population_template with rigid registration
  population_template \
    ${fba}/template/intra-temps/${sub}/fods/ \
    ${fba}/template/intra-temps/${sub}/${sub}_avg.mif \
    -mask_dir ${fba}/template/intra-temps/${sub}/masks/ \
    -voxel_size 1.25 \
    -type rigid \
    -linear_transformations_dir ${fba}/template/intra-temps/${sub}/xfms

  echo "Done with subject: ${sub}"
  echo "----------------------------"

done

# Now generate the study template
template=${der}/study_files/fba/template.txt
while read -r sub; do
  mkdir -p ${fba}/template/study_template/fods
  ln -sf ${fba}/template/intra-temps/${sub}/${sub}_avg.mif ${fba}/template/study_template/fods/
done < $template

population_template \
  ${fba}/template/study_template/fods/ \
  ${fba}/template/study_template/wmfod_template.mif \
  -voxel_size 1.25
  
# Population_template by default is: rigid_affine_nonlinear
# Want to use this on all subjects for consistency

mkdir -p \
${fba}/template/study_template/intra-template_xfms

for dir in ${fba}/data/sub-*; do
  sub=$(basename ${dir})

  mrregister \
  ${fba}/template/intra-temps/${sub}/${sub}_avg.mif \
  -type rigid_affine_nonlinear \
  ${fba}/template/study_template/wmfod_template.mif \
  -nl_warp \
  ${fba}/template/study_template/intra-template_xfms/${sub}_intra-temp_warp.mif \
  ${fba}/template/study_template/intra-template_xfms/tmp_${sub}_warp.mif # remove these

  rm ${fba}/template/study_template/intra-template_xfms/tmp_${sub}_warp.mif

  # Now register and transform all native timepoints
  # Use the affine to intra-pop and the warp of intra-pop to group

  for ses in ses-01 ses-02 ses-03; do

  # Compose for ease and later when need single file
  transformcompose \
  ${fba}/template/intra-temps/${sub}/xfms/${sub}_${ses}.txt \
  ${fba}/template/study_template/intra-template_xfms/${sub}_intra-temp_warp.mif \
  ${dir}/${ses}/fod/${sub}-template_warp.mif

  # Apply to the mask
  mrtransform \
    ${dir}/${ses}/fod/${sub}_b0_brain_mask_us.mif \
    -warp ${dir}/${ses}/fod/${sub}-template_warp.mif \
    -interp nearest -datatype bit \
    ${dir}/${ses}/fod/${sub}_b0_mask_us-template.mif \
    -template ${fba}/template/study_template/wmfod_template.mif

  done
done

# Compute wmfod template mask
mrmath \
${fba}/data/*/*/fod/*_b0_mask_us-template.mif \
min \
${fba}/template/study_template/template_mask.mif \
-datatype bit

# Define group white matter fixel mask and estimate fixel metrics for each patient
fod2fixel
-mask ${fba}/template/study_template/template_mask.mif \
-fmls_peak_value 0.08 \
${fba}/template/study_template/wmfod_template.mif \
${fba}/template/study_template/fixel_mask_08

# Segment FOD images to estimate fixels and their AFD
mkdir ${fba}/template/study_template/metrics
for dir in ${fba}/data/sub-*; do
  sub=$(basename ${dir})
  for ses in ses-01 ses-02 ses-03; do

    mkdir -p ${dir}/${ses}/fixels

    mrtransform ${dir}/${ses}/fod/${sub}_wmfod.mif \
    -warp ${dir}/${ses}/fod/${sub}-template_warp.mif \
    --reorient_fod no \
    ${dir}/${ses}/fod/${sub}_wmfod_noreo.mif
    
    fod2fixel \
    -mask ${fba}/template/study_template/template_mask.mif \
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
    ${dir}/${ses}/fixels/${sub}_fixel-template/${ses}_fd.mif \
    ${fba}/template/study_template/fixel_mask_08 \
    ${fba}/template/study_template/metrics/fd \
    ${sub}_${ses}.mif

    # Compute FC metric
    warp2metric \
    ${dir}/${ses}/fod/${sub}-template_warp.mif \
    -fc ${fba}/template/study_template/fixel_mask_08 \
    ${fba}/template/study_template/metrics/fc \
    ${sub}_${ses}.mif \
    -force
  done
done

# Copy files for, and compute log-fc
mkdir ${fba}/template/study_template/metrics/log_fc
cp ${fba}/template/study_template/metrics/fc/index.mif ${fba}/template/study_template/metrics/fc/directions.mif \
${fba}/template/study_template/metrics/log_fc

for dir in ${fba}/data/sub-*; do
   sub=$(basename ${dir})
   for ses in ses-01 ses-02 ses-03; do
     mrcalc \
     ${fba}/template/study_template/metrics/fc/${sub}_${ses}.mif \
     -log ${fba}/template/study_template/metrics/log_fc/${sub}_${ses}.mif \
     -force
   done
done

# Copy files for, and compute fdc
mkdir ${fba}/template/study_template/metrics/fdc
cp ${fba}/template/study_template/metrics/fc/index.mif \
${fba}/template/study_template/metrics/fc/directions.mif \
${fba}/template/study_template/metrics/fdc

for dir in ${fba}/data/sub-*; do
   sub=$(basename ${dir})
   for ses in ses-01 ses-02 ses-03; do
     mrcalc \
     ${fba}/template/study_template/metrics/fd/${sub}_${ses}.mif \
     ${fba}/template/study_template/metrics/fc/${sub}_${ses}.mif -mult \
     ${fba}/template/study_template/metrics/fdc/${sub}_${ses}.mif \
     -force
   done
done

# Generate tractogram
tckgen \
-angle 22.5 \
-maxlen 250 \
-minlen 10 \
-power 1.0 \
${fba}/template/study_template/wmfod_template.mif \
-seed_image ${fba}/template/study_template/template_mask.mif \
-select 20000000 \
-cutoff 0.08  \ # matched to the fod2fixel threshold
${fba}/template/study_template/tractogram_20mil.tck

# SIFT tractogram
tcksift \
${fba}/template/study_template/tractogram_20mil.tck \
${fba}/template/study_template/wmfod_template.mif \
${fba}/template/study_template/tractogram_2mil_SIFT.tck \
-term_number 2000000

# Compute fixel-fixel connectivity matrix
fixelconnectivity \
${fba}/template/study_template/fixel_mask_08/ \
${fba}/template/study_template/tractogram_2mil_SIFT.tck \
${fba}/template/study_template/matrix/ \
-force

# Smooth metric data using the fixel-fixel connectivity matrix
for metric in fd log_fc fdc; do
fixelfilter \
${fba}/template/study_template/metrics/${metric} \
smooth ${fba}/template/study_template/metrics/${metric}_smooth \
-matrix ${fba}/template/study_template/matrix/
done

# Generate wmfod<->MNI (0.5mm) xfm
mrconvert \
${fba}/template/study_template/wmfod_template.mif \
-coord 3 0 -axes 0,1,2 \
${fba}/template/study_template/wmfod_template.nii.gz

mri_synthstrip \
-i ${fba}/template/study_template/wmfod_template.nii.gz \
-o ${fba}/template/study_template/wmfod_template_brain.nii.gz \
-m ${fba}/template/study_template/wmfod_template_brain_mask.nii.gz

MNI_T2=/Applications/leaddbs/templates/space/MNI152NLin2009bAsym/t2.nii
antsRegistrationSyN.sh \
  -d 3 \
  -f ${MNI_T2} \
  -m ${fba}/template/study_template/wmfod_template.nii.gz \
  -o ${fba}/template/study_template/wmfod_template-MNI_ \
  -n 12

# Group analyses:

# Percent change images for postop. timepoints
fdcdir=${fba}/template/study_template/metrics/fdc_smooth
analysis=${fba}/analysis
# Create output directories
mkdir -p ${analysis}/group_fba/6m_pc \
${analysis}/group_fba/12m_pc

# Loop through ses-01 files
for pre in ${fdcdir}/sub-*_ses-01.mif; do
  # Extract subject ID (e.g., sub-001)
  sub=$(basename $pre | cut -d '_' -f 1)

  # Define paths to sessions
  fdc_6m=${fdcdir}/${sub}_ses-02.mif
  fdc_12m=${fdcdir}/${sub}_ses-03.mif

  # Output filenames
  pc_6m=${analysis}/group_fba/6m_pc/${sub}.mif
  pc_12m=${analysis}/group_fba/12m_pc/${sub}.mif

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

# 6-month and 12-month percent change analyses
for tp in 6m 12m; do
  ln -s ${fba}/template/study_template/metrics/fdc_smooth/directions.mif \
  ${fba}/template/study_template/metrics/fdc_smooth/index.mif \
  ${analysis}/group_fba/${tp}_pc/

  fixelcfestats \
  ${analysis}/group_fba/${tp}_pc \
  ${analysis}/cfe_files/group_analyses/${tp}_subs.txt \
  ${analysis}/cfe_files/group_analyses/${tp}_demeaned.txt \
  ${analysis}/cfe_files/group_analyses/corr_con.txt \
  ${fba}/template/study_template/matrix \
  ${analysis}/cfe_files/group_analyses/${tp}_results
done

# Preoperative severity analysis
fixelcfestats \
  ${fba}/template/study_template/metrics/fdc_smooth  \
  ${analysis}/cfe_files/group_analyses/pre_subs.txt \
  ${analysis}/cfe_files/group_analyses/pre_demeaned.txt \
  ${analysis}/cfe_files/group_analyses/corr_con.txt \
  ${fba}/template/study_template/matrix \
  ${analysis}/cfe_files/group_analyses/pre_results

###################

# Store the FDC value of the sweetspot per session (else NA)

# BUT NEED TO FIRST FLIP THE SWEETSPOT ONTO OTHER HEMI AND THEN WORK OUT WAY OF INDEXING THE TARGETED HEMI TO MULTIPLY BY WHICH HEMI SWEETSPOT

echo "sub-id,ses-01,ses-02,ses-03" > ${fba}/analysis/sweetspot/lh.csv 
echo "sub-id,ses-01,ses-02,ses-03" > ${fba}/analysis/sweetspot/rh.csv
sweetspot=${fba}/analysis/sweetspot/sweetspot.nii.gz

# Loop through each subject
for dir in ${fba}/data/sub-*; do
  sub=$(basename ${dir})
  echo "Processing $sub"

  lh_line="$sub"
  rh_line="$sub"

  # Loop through sessions
  for ses in ses-01 ses-02 ses-03; do
    fdc_mif="${fba}/template/study_template/metrics/fdc_smooth/${sub}_${ses}.mif"
    fdc_nii="${fba}/data/${sub}/${ses}/fod/${sub}_fdc_map.nii.gz"

    # Check if input file exists
    if [[ -f "$fdc_mif" ]]; then
      fixel2voxel "$fdc_mif" mean "$fdc_nii"

      for hemi in lh rh; do
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

# Hypointensity permutation modelling
for dir in ${fba}/data/sub-*; do
  sub=$(basename ${dir})

  for hemi in lh rh; do

  # 6-months
  fslmaths \
    ${fba}/data/${sub}/ses-01/fod/${sub}_fdc_${hemi}_hypo.nii.gz \
    -sub ${fba}/data/${sub}/ses-02/fod/${sub}_fdc_${hemi}_hypo.nii.gz \
    -div ${fba}/data/${sub}/ses-01/fod/${sub}_fdc_${hemi}_hypo.nii.gz \
    -mul 100 ${fba}/data/${sub}/ses-02/${sub}_hypo_pc.nii.gz

  # 12-months
  fslmaths \
    ${fba}/data/${sub}/ses-01/fod/${sub}_fdc_${hemi}_hypo.nii.gz \
    -sub ${fba}/data/${sub}/ses-03/fod/${sub}_fdc_${hemi}_hypo.nii.gz \
    -div ${fba}/data/${sub}/ses-01/fod/${sub}_fdc_${hemi}_hypo.nii.gz \
    -mul 100 ${fba}/data/${sub}/ses-03/${sub}_hypo_pc.nii.gz

  done
done

fslrandomise...








# Calculate brain atrophy using a regression-residuals method (TBV~ICV)


