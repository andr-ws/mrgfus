#! /bin/bash

export FREESURFER_HOME=/Applications/freesurfer
source $FREESURFER_HOME/SetUpFreeSurfer.sh

# Compute FOD response functions for DWI
base=~/imaging/datasets/mrgfus
raw=${base}/rawdata
dwi=${base}/derivatives/data/dwi
fba=${base}/derivatives/fba
sf=${base}/derivatives/study_files

hemis=${sf}/fba/hemis.txt # sub-id hemi
while read -r sub hemi; do

  for ses in ses-01 ses-02 ses-03; do
    eddy=${dwi}/${sub}/${ses}/eddy
    fod=${fba}/data/${sub}/${ses}/fod
    
    if [ ! -d ${eddy} ]; then
      echo "Skipping subject ${sub} ${ses} - this directory is not found!"
      continue
    fi
  
    echo "Processing subject ${sub} for session ${ses}"
    mkdir -p ${fod}

    # Need to check if here if a subject is rhetorical targeted
    if [[ "$hemi" == "lh" ]]; then
    
      # Flip x-axis
      mrtransform \
        ${eddy}/${sub}_dwi_edc.nii.gz \
        ${fod}/tmp_${sub}_dwi.mif \
        -flip 0 \
        -fslgrad \
        ${eddy}/${sub}_dwi_edc.eddy_rotated_bvecs \
        ${dwi}/${sub}/${ses}/${sub}_${ses}_acq-dwi.bval
    elif [[ "$hemi" == "rh" ]]; then
    
      mrconvert \
        ${eddy}/${sub}_dwi_edc.nii.gz \
        ${fod}/tmp_${sub}_dwi.mif \
        -fslgrad \
        ${eddy}/${sub}_dwi_edc.eddy_rotated_bvecs \
        ${dwi}/${sub}/${ses}/${sub}_${ses}_acq-dwi.bval
    else
      echo "ERROR: unexpected value for hemi: '$hemi'" >&2
      exit 1
    fi
    
    # Bias correct
    dwibiascorrect \
      ants \
      ${fod}/tmp_${sub}_dwi.mif \
      ${fod}/fod/tmp_${sub}_dwi_dn.mif

    # Compute average tissue-response function
    dwi2response \
      dhollander \
      ${fod}/fod/tmp_${sub}_dwi_dn.mif \
      ${fod}/fod/response_wm.txt \
      ${fod}/${ses}/fod/response_gm.txt \
      ${fod}/fod/response_csf.txt

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
      -m ${fod}/tmp_${sub}_b0_brain_mask_us.nii.gz

    # Convert brain mask to mif format
    mrconvert \
      ${fod}/tmp_${sub}_b0_brain_mask_us.nii.gz \
      ${fod}/${sub}_b0_brain_mask_us.mif

    # Clean up
    rm ${fod}/tmp*
      
  done
done < ${hemis}

# Create group tissue response functions
mkdir -p ${fba}/rfuncs

for tissue in wm gm csf; do
  responsemean \
    ${fba}/*/*/*/fod/response_${tissue}.txt \
    ${fba}/rfuncs/group_response_${tissue}.txt
done

# Create individual tissue response functions
for dir in ${fba}/data/sub-*; do
  sub=$(basename ${dir})

  for ses in ses-01 ses-02 ses-03; do
    fod=${dir}/${ses}/fod
    
    ss3t_csd_beta1 \
      ${fod}/${sub}_dwi_us.mif \
      ${fba}/rfuncs/group_response_wm.txt \
      ${fod}/tmp_${sub}_wmfod.mif \
      ${fba}/rfuncs/group_response_gm.txt \
      ${fod}/tmp_${sub}_gmfod.mif \
      ${fba}/rfuncs/group_response_csf.txt \
      ${fod}/tmp_${sub}_csffod.mif \
      -mask ${fod}/${sub}_b0_brain_mask_us.mif

    mtnormalise \
      ${fod}/tmp_${sub}_wmfod.mif ${fod}/${sub}_wmfod.mif \
      ${fod}/tmp_${sub}_gmfod.mif ${fod}/${sub}_gmfod.mif \
      ${fod}/tmp_${sub}_csffod.mif ${fod}/${sub}_csffod.mif \
      -mask ${fod}/${sub}_b0_brain_mask_us.mif

    # Clean up
    rm \
      ${fod}/${sub}_dwi_us.mif \
      ${fod}/tmp_${sub}_wmfod.mif \
      ${fod}/tmp_${sub}_gmfod.mif \
      ${fod}/tmp_${sub}_csffod.mif

  done
done

# Longitudinal FBA (oh god)!

# Desc: generate rigid intra-pop templates and affine xfms for each session
# Generate a study template and write out nlin warps for those individuals
# Generate nonlinear study template warps for those not in the template
# Compose the linear and nonlinear transforms

for dir in ${fba}/data/sub-*; do
  sub=$(basename ${dir})

  itemp=${fba}/template/itemps/${sub}
  mkdir -p ${itemp}/fods ${itemp}/masks

  # Loop through possible sessions
  for ses in ses-01 ses-02 ses-03; do
    fod=${fba}/data/${sub}/${ses}/fod
    
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
  population_template \
    ${itemp}/fods/ \
    ${itemp}/fods/${sub}_itemp.mif \
    -mask_dir ${itemp}/masks/ \
    -voxel_size 1.25 \
    -type rigid \
    -linear_transformations_dir ${itemp}/xfms

  # Rename for clarity
  mv ${itemp}/xfms/${sub}_${ses}.txt ${itemp}/xfms/${sub}_${ses}-itemp_rigid.txt 

  echo "Done with subject: ${sub}"
  echo "----------------------------"

done

# Now generate the study template
template=${sf}/fba/template.txt
mkdir -p ${fba}/template/study_template/fods

while read -r sub; do
  ln -s \
    ${fba}/template/itemps/${sub}/fods/${sub}_itemp.mif \
    ${fba}/template/study_template/fods/
done < $template

wmfod=${fba}/template/study_template/wmfod
mkdir -p ${wmfod}

population_template \
  ${fba}/template/study_template/fods/ \
  ${wmfod}/wmfod_template.mif \
  -voxel_size 1.25

for dir in ${fba}/data/sub-*; do
  sub=$(basename ${dir})

  itemp=${fba}/template/itemps
  
  mrregister \
    ${itemp}/${sub}/fods/${sub}_itemp.mif \
    -type rigid_affine_nonlinear \
    ${wmfod}/wmfod_template.mif \
    -nl_warp \ # rigid_affine_nonlinear
    ${itemp}/${sub}/xfms/${sub}_itemp-temp_warp.mif \
    ${itemp}/${sub}/xfms/tmp_${sub}.mif # remove these

  rm ${itemp}/${sub}/xfms/tmp_${sub}.mif

  # Now register and transform all native timepoints
  # Use the affine to intra-pop and the warp of intra-pop to group

  for ses in ses-01 ses-02 ses-03; do

    fod=${dir}/${ses}/fod

    # Compose for ease and later when need single file
    transformcompose \
      ${itemp}/${sub}/xfms/${sub}_${ses}_itemp_rigid.txt \
      ${itemp}/${sub}/xfms/${sub}_itemp-temp_warp.mif \
      ${fod}/${sub}-template_warp.mif

    # Apply to the mask
    mrtransform \
      ${fod}/${sub}_b0_brain_mask_us.mif \
      -warp ${fod}/${sub}-template_warp.mif \
      -interp nearest -datatype bit \
      ${fod}/${sub}_b0_mask_us-template.mif \
      -template ${wmfod}/wmfod_template.mif
  done
done

# Compute wmfod template mask
mrmath \
  ${fba}/data/*/*/fod/*_b0_mask_us-template.mif \
  min \
  ${wmfod}/template_mask.mif \
  -datatype bit

# Define group white matter fixel mask and estimate fixel metrics for each patient
fod2fixel
  -mask ${wmfod}/template_mask.mif \
  -fmls_peak_value 0.08 \
  ${wmfod}/wmfod_template.mif \
  ${fba}/template/study_template/fixel_mask_08

# Segment FOD images to estimate fixels and their AFD
unsmoothed=${fba}/template/study_template/metrics/unsmoothed
mkdir -p ${unsmoothed}

for dir in ${fba}/data/sub-*; do
  sub=$(basename ${dir})
  for ses in ses-01 ses-02 ses-03; do

    fod=${dir}/${ses}/fod
    fix=${dir}/${ses}/fixels

    mkdir -p ${fix}

    mrtransform \
      ${fod}/${sub}_wmfod.mif \
      -warp ${fod}/${sub}-template_warp.mif \
      --reorient_fod no \
      ${fod}/${sub}_wmfod_noreo.mif
    
    fod2fixel \
      -mask ${wmfod}/template_mask.mif \
      ${fod}/${sub}_wmfod_noreo.mif \
      ${fix}/fixel-template_noreo \
      -afd ${sub}_fd.mif
   
    # Reorient fixels
    fixelreorient \
      ${fix}/fixel-template_noreo \
      ${fod}/${sub}-template_warp.mif \
      ${dir}/${ses}/fixels/fixel-template

    # Assign subjects fixels to template fixels
    fixelcorrespondence \
      ${fix}/fixel-template/${sub}_fd.mif \
      ${fba}/template/study_template/fixel_mask_08 \
      ${unsmoothed}/fd \
      ${sub}_${ses}.mif

    # Compute FC metric
    warp2metric \
      ${fod}/${sub}-template_warp.mif \
      -fc ${fba}/template/study_template/fixel_mask_08 \
      ${unsmoothed}/fc \
      ${sub}_${ses}.mif
  done
done

# Copy files for, and compute log-fc
mkdir ${unsmoothed}/log_fc
ln -s \
  ${unsmoothed}/fc/index.mif \
  ${unsmoothed}/fc/directions.mif \
  ${unsmoothed}/log_fc

for dir in ${fba}/data/sub-*; do
   sub=$(basename ${dir})
   for ses in ses-01 ses-02 ses-03; do
     mrcalc \
       ${unsmoothed}/fc/${sub}_${ses}.mif \
       -log ${unsmoothed}/log_fc/${sub}_${ses}.mif
   done
done

# Copy files for, and compute fdc
mkdir ${unsmoothed}/fdc
ln -s \
  ${unsmoothed}/fc/index.mif \
  ${unsmoothed}/fc/directions.mif \
  ${unsmoothed}/fdc

for dir in ${fba}/data/sub-*; do
   sub=$(basename ${dir})
   for ses in ses-01 ses-02 ses-03; do
     mrcalc \
       ${unsmoothed}/fd/${sub}_${ses}.mif \
       ${unsmoothed}/fc/${sub}_${ses}.mif -mult \
       ${unsmoothed}/fdc/${sub}_${ses}.mif
   done
done

tracts=${fba}/template/study_template/tracts
mkdir -p ${tracts}

# Generate tractogram
tckgen \
  -angle 22.5 \
  -maxlen 250 \
  -minlen 10 \
  -power 1.0 \
  ${wmfod}/wmfod_template.mif \
  -seed_image ${wmfod}/template_mask.mif \
  -select 20000000 \
  -cutoff 0.08  \ # matched to the fod2fixel threshold
  ${tracts}/tracts-20m.tck

# SIFT tractogram
tcksift \
  ${tracts}/tracts-20m.tck \
  ${wmfod}/wmfod_template.mif \
  ${tracts}/tracts-2m.tck \
  -term_number 2000000

# Compute fixel-fixel connectivity matrix
fixelconnectivity \
  ${fba}/template/study_template/fixel_mask_08/ \
  ${tracts}/tracts-2m.tck \
  ${fba}/template/study_template/matrix/

# Smooth metric data using the fixel-fixel connectivity matrix
smoothed=${fba}/template/study_template/metrics/smoothed
for metric in fd log_fc fdc; do
  mkdir -p ${smoothed}/${metric}
  
  fixelfilter \
    ${unsmoothed}/${metric} \
    smooth ${smoothed}/${metric} \
    -matrix ${fba}/template/study_template/matrix/
done

# Generate wmfod<->MNI (0.5mm) xfm
# USE SYNTHMORPH HERE>>
#mrconvert \
#${fba}/template/study_template/wmfod_template.mif \
#-coord 3 0 -axes 0,1,2 \
#${fba}/template/study_template/wmfod_template.nii.gz

#mri_synthstrip \
#-i ${fba}/template/study_template/wmfod_template.nii.gz \
#-o ${fba}/template/study_template/wmfod_template_brain.nii.gz \
#-m ${fba}/template/study_template/wmfod_template_brain_mask.nii.gz

#MNI_T2=/Applications/leaddbs/templates/space/MNI152NLin2009bAsym/t2.nii
#antsRegistrationSyN.sh \
#  -f ${MNI_T2} \
#  -d 3 \
#  -m ${fba}/template/study_template/wmfod_template.nii.gz \
#  -n 12
#  -o ${fba}/template/study_template/wmfod_template-MNI_ \

# Group analyses:

# Percent change images for postop. timepoints
fdcdir=${fba}/template/study_template/metrics/smoothed/fdc
analysis=${fba}/analysis
# Create output directories
mkdir -p ${analysis}/group_fba/6m_pc \
${analysis}/group_fba/12m_pc \
${analysis}/group_fba/pre

# Loop through ses-01 files
for pre in ${fdcdir}/sub-*_ses-01.mif; do
  # Extract subject ID (e.g., sub-001)
  sub=$(basename $pre | cut -d '_' -f 1)

  # Populate preop directory
  fdc_pre=${fdcdir}/${sub}_ses-01.mif
  ln -s ${fdc_pre} ${analysis}/group_fba/pre/${sub}.mif

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

# Preoperative severity analysis
ln -s ${fba}/template/study_template/metrics/smoothed/fdc/directions.mif \
  ${fba}/template/study_template/metrics/smoothed/fdc/index.mif \
  ${analysis}/group_fba/pre/
  
fixelcfestats \
  ${analysis}/group_fba/pre \
  ${analysis}/cfe_files/group_analyses/pre_subs.txt \
  ${analysis}/cfe_files/group_analyses/pre_demeaned.txt \
  ${analysis}/cfe_files/group_analyses/corr_con.txt \
  ${fba}/template/study_template/matrix \
  ${analysis}/cfe_files/group_analyses/pre_results

# 6-month and 12-month percent change analyses
for tp in 6m 12m; do
  ln -s ${fba}/template/study_template/metrics/smoothed/fdc/directions.mif \
  ${fba}/template/study_template/metrics/smoothed/fdc/index.mif \
  ${analysis}/group_fba/${tp}_pc/

  fixelcfestats \
  ${analysis}/group_fba/${tp}_pc \
  ${analysis}/cfe_files/group_analyses/${tp}_subs.txt \
  ${analysis}/cfe_files/group_analyses/${tp}_demeaned.txt \
  ${analysis}/cfe_files/group_analyses/corr_con.txt \
  ${fba}/template/study_template/matrix \
  ${analysis}/cfe_files/group_analyses/${tp}_results
done

###################

# Store the FDC value of the sweetspot per session (else NA)
mean_fdc_csv=${fba}/analysis/sweetspot/sweetspot_mean-fdc.csv
echo "sub-id,ses-01,ses-02,ses-03" > ${mean_fdc_csv}
mkdir -p ${fba}/analysis/sweetspot/mean_fdc_maps

# Loop through each subject
while read -r sub; do
  echo "Processing $sub"
  sweetspot=POINT_2_WMFOD_SS

  # Initialize CSV line
  line="$sub"

  # Loop through sessions
  for ses in ses-01 ses-02 ses-03; do
    fdc_mif=${fba}/template/study_template/metrics/smoothed/fdc/${sub}_${ses}.mif
    fdc_map=${fba}/analysis/sweetspot/mean_fdc_maps/${sub}_${ses}.nii.gz

    # Check if input file exists
    if [[ -f "$fdc_mif" ]]; then
      fixel2voxel \
        ${fdc_mif} \
        mean \
        ${fdc_map} \
        -weighted <UNSURE OF WHAT THIS NEED TO BE> ###########################
        
      fslmaths \
        ${fdc_map} \
        -mul \
        ${sweetspot} \
        ${fdc_map} # overwrite
        
      value=$(fslstats "$fdc_map" -M)
    else
      echo "  $ses: Missing FDC file."
      value="NA"
    fi

    # Append value to line
    line="${line},${value}"
  done
  # Append full row to CSV
  echo ${line} >> ${mean_fdc_csv}
done < ${hemis}

# Calculate brain atrophy using a regression-residuals method (TBV~ICV)
