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
    if [[ "$hemi" == "rh" ]]; then
    
      # Flip x-axis
      mrtransform \
        ${eddy}/${sub}_dwi_edc.nii.gz \
        ${fod}/tmp_${sub}_dwi.mif \
        -flip 0 \
        -fslgrad \
        ${eddy}/${sub}_dwi_edc.eddy_rotated_bvecs \
        ${dwi}/${sub}/${ses}/${sub}_${ses}_acq-dwi.bval
        
    elif [[ "$hemi" == "lh" ]]; then
    
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
      ${fod}/tmp_${sub}_dwi_us.mif

    # Extract upsampled b0
    dwiextract \
      ${fod}/tmp_${sub}_dwi_us.mif \
      - -bzero | mrmath - mean ${fod}/tmp_${sub}_b0_us.mif \
      -axis 3

    mrconvert \
      ${fod}/tmp_${sub}_b0_us.mif \
      ${fod}/tmp_${sub}_b0_us.nii.gz

    # Brain extract and create mask
    # Keep the us_b0 for later (MNI to template warping)
    mri_synthstrip \
      -i ${fod}/tmp_${sub}_b0_us.nii.gz \
      -o ${fod}/${sub}_b0_brain_us.nii.gz \
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

    source activate ~/mrtrix_env/bin/activate
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

  itemp=${fba}/template/intra-temps/${sub}
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
  # Was -type rigid , now running rigid_affine_nonlinear, as default
  population_template \
    ${itemp}/fods/ \
    ${itemp}/fods/${sub}_itemp.mif \
    -mask_dir ${itemp}/masks/ \
    -voxel_size 1.25
    #-warp_dir ${itemp}/xfms - cancelled this as doesn't recognise as
    # warp fields. These are recomputed later.

  mrconvert \
    ${itemp}/fods/${sub}_itemp.mif \
    -coord 3 0 -axes 0,1,2 \
    ${itemp}/fods/tmp_${sub}_itemp.nii.gz

  fslmaths \
    ${itemp}/fods/tmp_${sub}_itemp.nii.gz \
    -nan \
    ${itemp}/fods/tmp_${sub}_itemp_no-nan.nii.gz

  mri_synthstrip \
    -i ${itemp}/fods/tmp_${sub}_itemp_no-nan.nii.gz \
    -o ${itemp}/fods/tmp_${sub}_itemp_strip.nii.gz \
    -m ${itemp}/fods/tmp_${sub}_itemp_mask.nii.gz \
    -t 12

  mrconvert \
    ${itemp}/fods/tmp_${sub}_itemp_mask.nii.gz \
    ${itemp}/fods/${sub}_itemp_mask.mif

  rm ${itemp}/fods/tmp*

  # Rename for clarity
  #mv ${itemp}/xfms/${sub}_${ses}.txt ${itemp}/xfms/${sub}_${ses}-itemp_rigid.txt

  echo "Done with subject: ${sub}"
  echo "----------------------------"

done

# Now generate the study template
template=${sf}/fba/template.txt
mkdir -p ${fba}/template/study_template/fods ${fba}/template/study_template/masks

while read -r sub; do
  ln -s ${fba}/template/intra-temps/${sub}/fods/${sub}_itemp.mif ${fba}/template/study_template/fods/
  ln -s ${fba}/template/intra-temps/${sub}/fods/${sub}_itemp_mask.mif ${fba}/template/study_template/masks/${sub}_itemp.mif
done < $template

wmfod=${fba}/template/study_template/wmfod
mkdir -p ${wmfod}

population_template \
  ${fba}/template/study_template/fods/ \
  ${wmfod}/wmfod_template.mif \
  -mask_dir ${fba}/template/study_template/masks/ \
  -voxel_size 1.25

for dir in ${fba}/data/sub-*; do
  sub=$(basename ${dir})

  itemp=${fba}/template/intra-temps

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
    # Perhaps need to ammend to compose two warps
    
    transformcompose \
      ${itemp}/${sub}/xfms/${sub}_${ses}.mif \ # The warp from native to intra
      ${itemp}/${sub}/xfms/${sub}_itemp-temp_warp.mif \ # The warp from intra to template
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

# Generate wmfod to MNI (0.5mm) xfm - T2 the best after testing
T2=/Applications/leaddbs/templates/space/MNI152NLin2009bAsym/t2.nii.gz
mrconvert \
${fba}/template/study_template/wmfod_template.mif \
  -coord 3 0 -axes 0,1,2 \
  ${fba}/template/study_template/wmfod_template.nii.gz

fslmaths \
  ${fba}/template/study_template/wmfod_template.nii.gz \
  -nan \
  ${fba}/template/study_template/wmfod_template_noNaN.nii.gz

mri_synthstrip \
  -i ${fba}/template/study_template/wmfod/wmfod_template_noNaN.nii.gz \
  -o ${fba}/template/study_template/wmfod/wmfod_template_noNaN_brain.nii.gz \
  -m ${fba}/template/study_template/wmfod/wmfod_template_noNaN_brain_mask.nii.gz

mri_synthmorph \
  register \
  -o ${fba}/template/study_template/wmfod/wmfod_template-MNI.nii.gz \
  -O ${fba}/template/study_template/wmfod/MNI-wmfod_template.nii.gz \
  ${fba}/template/study_template/wmfod/wmfod_template_noNaN.nii.gz \
  ${T2} \
  -t ${fba}/template/study_template/wmfod/wmfod_template-MNI_warp.nii.gz \
  -T ${fba}/template/study_template/wmfod/MNI-wmfod_template_warp.nii.gz

# Group analyses:
metricbase=${fba}/template/study_template/metrics/smoothed
analysis=${fba}/analysis

# Create output directories
# Loop through ses-01 files
for metric in fdc fd log_fc; do

  mkdir -p 
  ${analysis}/group_fba/${metric}/preop \
  ${analysis}/group_fba/${metric}/6m \
  ${analysis}/group_fba/${metric}/12m

  for pre in ${metricbase}/${metric}/sub-*_ses-01.mif; do
    metricdir=${metricbase}/${metric}
    
    # Extract subject ID (e.g., sub-001)
    sub=$(basename $pre | cut -d '_' -f 1)

    # Populate preop directory
    metric_pre=${metricdir}/${sub}_ses-01.mif
    ln -s ${metric_pre} ${analysis}/group_fba/preop/${sub}.mif

    # Define paths to sessions
    metric_6m=${metricdir}/${sub}_ses-02.mif
    metric_12m=${metricdir}/${sub}_ses-03.mif

    metric_diff_6m=${analysis}/group_fba/${metric}/6m/${sub}.mif
    metric_diff_12m=${analysis}/group_fba/${metric}/12m/${sub}.mif

    # Compute 6-month percent change if ses-02 exists
    if [ -f $fdc_6m ]; then
      echo "Computing 6m difference maps for ${sub}"
      mrcalc $pre $fdc_6m -subtract $diff_6m
    else
      echo "Skipping 6m for ${sub} (missing ses-02)"
    fi

    # Compute 12-month percent change if ses-03 exists
    if [ -f $fdc_12m ]; then
      echo "Computing 12m difference maps for ${sub}"
      mrcalc $pre $fdc_12m -subtract $diff_12m
    else
      echo "Skipping 12m for ${sub} (missing ses-03)"
    fi
  done
done




analysis=${fba}/analysis

# Preoperative severity analyses
for metric in fdc log_fc fd; do
  ln -s ${fba}/template/study_template/metrics/smoothed/${metric}/directions.mif \
    ${fba}/template/study_template/metrics/smoothed/${metric}/index.mif \
    ${fba}/template/study_template/metrics/preop/${metric}/
done

# Log_fc and FDC
for metric in log_fc fdc; do
  fixelcfestats \
    ${fba}/template/study_template/metrics/preop/${metric}/ \
    ${analysis}/cfe_files/group_analyses/preop/preop_subs.txt \
    ${analysis}/cfe_files/group_analyses/preop/unit_fc.txt \
    ${analysis}/cfe_files/group_analyses/contrasts/fc/preop.txt \
    ${fba}/template/study_template/matrix \
    ${analysis}/cfe_files/group_analyses/preop/${metric}
done

# FD
fixelcfestats \
  ${fba}/template/study_template/metrics/preop/fd \
  ${analysis}/cfe_files/group_analyses/preop/preop_subs.txt \
  ${analysis}/cfe_files/group_analyses/preop/unit_fd.txt \
  ${analysis}/cfe_files/group_analyses/contrasts/fd/preop.txt \
  ${fba}/template/study_template/matrix \
  ${analysis}/cfe_files/group_analyses/preop/fd


# 6-month and 12-month difference analyses
for metric in log_fc fdc fd; do
  for tp in 6m 12m; do
    ln -s ${fba}/template/study_template/metrics/smoothed/${metric}/directions.mif \
      ${fba}/template/study_template/metrics/smoothed/${metric}/index.mif \
      ${fba}/template/study_template/metrics/difference/${metric}/${tp}/
  done
done

# log_fc and FDC
for tp in 6m 12m; do
  for metric in log_fc fdc; do
    fixelcfestats \
      ${fba}/template/study_template/metrics/difference/${metric}/${tp}/ \
      ${analysis}/cfe_files/group_analyses/${tp}/${tp}_subs.txt \
      ${analysis}/cfe_files/group_analyses/${tp}/unit_fc.txt \
      ${analysis}/cfe_files/group_analyses/contrasts/fc/diff.txt \
      ${fba}/template/study_template/matrix \
      ${analysis}/cfe_files/group_analyses/${tp}/${metric}
  done

  # FD
  fixelcfestats \
    ${fba}/template/study_template/metrics/difference/${metric}/${tp}/ \
    ${analysis}/cfe_files/group_analyses/${tp}/${tp}_subs.txt \
    ${analysis}/cfe_files/group_analyses/${tp}/unit_fd.txt \
    ${analysis}/cfe_files/group_analyses/contrasts/fd/diff.txt \
    ${fba}/template/study_template/matrix \
    ${analysis}/cfe_files/group_analyses/${tp}/fd
done
