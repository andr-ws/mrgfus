#! /bin/bash

# DESC: This code performs registration of native lesion T2s to MNI space,
# non-linearly flips those on the right hemisphere (to the left side), 
# performs statistical mapping and also warps the sweetspot into wmfod space.

# Note, this requires all lesions to already be manually defined...

raw=/Volumes/LA_4TB/datasets/mrgfus/rawdata
der=/Volumes/LA_4TB/datasets/mrgfus/derivatives
sf=${der}/study_files
lesions=${der}/lesions
MNI=/Applications/leaddbs/templates/space/MNI152NLin2009bAsym/t2.nii.gz

# Need to generate coreg btwn lesion T2 and ses-01 T1
for dir in ${raw}/sub-*; do
  sub=$(basename ${dir})
  
  mkdir -p ${lesions}/data/${sub}
  preT2=${dir}/rawdata/ses-01/anat/${sub}_acq-T2w.nii.gz
  immT2=${dir}/imm_postop/${sub}_acq-T2w.nii.gz
  preT1=${der}/anat/${sub}/ses-01/${sub}_ses-01_acq-T1w_biasco.nii.gz

  # Rigid transform
  antsRegistrationSyNQuick.sh \
    -d 3 \
    -f ${preT1} \
    -m ${preT2} \
    -o ${lesions}/data/${sub}/${sub}_T2w-ses-01_T1w_ \
    -t r \
    -n 12

  coregT2=${lesions}/data/${sub}/${sub}_T2w-ses-01_T1w_Warped.nii.gz

  # Register the immT2 to T2w in T1w space
  mri_synthmorph \
    register \
    ${immT2} \
    ${coregT2} \
    -o ${lesions}/data/${sub}/${sub}_immT2w-T1w.nii.gz \
    -t ${lesions}/data/${sub}/${sub}_immT2w-ses-01_T1w_1Warp.nii.gz \
    -T ${lesions}/data/${sub}/${sub}_immT2w-ses-01_T1w_1InverseWarp.nii.gz

  # Alternatively... imm-T2w straight to MNI T2!
  mri_synthmorph \
    register \
    ${immT2} \
    ${MNI}  \
    -o ${lesions}/data/${sub}/${sub}_immT2w-MNI.nii.gz \
    -t ${lesions}/data/${sub}/${sub}_immT2w-MNI_1Warp.nii.gz \
    -T ${lesions}/data/${sub}/${sub}_immT2w-MNI_1InverseWarp.nii.gz

  # Can then draw directly on MNI, or draw on original T2w and apply the transform
done

# Generate masks in native space
# ${lesions}/masks/nii/${sub}_lesion.nii.gz ...

for dir in ${lesions}/data/sub-*; do
  sub=$(basename ${dir})

  mkdir -p ${lesions}/masks/mni/${sub}

  # Apply the SynthMorph immT2-T1w and ANTs T1w-MNI warps
  antsApplyTransforms \
    -d 3 \
    -i ${lesions}/masks/nii/${sub}_lesion.nii.gz \
    -r ${MNI} \
    -t ${der}/anat/${sub}/ses-01/${sub}_T1w-MNI_05mm_1Warp.nii.gz \
    -t ${der}/anat/${sub}/ses-01/${sub}_T1w-MNI_05mm_0GenericAffine.mat \
    -t ${lesions}/data/${sub}/${sub}_immT2w-ses-01_T1w_1Warp.nii.gz \
    -o ${lesions}/masks/mni/${sub}/${sub}_lesion.nii.gz \
    --interpolation NearestNeighbor
done

# Flip rh lesions to lh
matlab \
  -nodisplay \
  -nosplash \
  -r "run('/Users/neuro-239/scripts/lesions_flip.m'); exit;"

# Modelling
mkdir ${lesions}/model
hemis=${sf}/fba/hemis.txt

mkdir ${lesions}/model/tmp
while read -r sub hemi; do
  if [[ $hemi == lh ]]; then
    ln -s \
    ${lesions}/masks/mni/${sub}/${sub}_lesion.nii.gz \
    ${lesions}/model/tmp/
  elif [[ $hemi == rh ]]; then
    ln -s \
    ${lesions}/masks/mni/${sub}/${sub}_flesion.nii.gz \
    ${lesions}/model/tmp/${sub}_lesion.nii.gz
  else 
    echo "Missing hemisphere information for ${sub}"
  fi
done < ${hemis}

# Generate files
# Create a 4D lesion input
mrcat \
  ${lesions}/model/tmp/*_lesion.nii.gz \
  ${lesions}/model/4d_lesions.nii.gz
# Create an N-map \
mrmath \
  ${lesions}/model/tmp/*lesion.nii.gz \
  sum \
  ${lesions}/model/lesion_n-map.nii.gz

# Tidy up tmp directory
rm -r ${lesions}/model/tmp

# Create study masks (remove lowest 10%, binarise and NaN 0)
fslmaths \
  ${lesions}/model/lesion_n-map.nii.gz \
  -thrP 10 \
  -bin \
  ${lesions}/model/lesion_n-map_thr.nii.gz

mrconvert \
  ${lesions}/model/lesion_n-map_thr.nii.gz \
  - | mrcalc - 1 nan -if \
  ${lesions}/model/lesion_n-map_thr.nii.g

# Perform sweetspot analysis
randomise \
  -i ${lesions}/model/4d_lesions.nii.gz \
  -o ${lesions}/model/randomise \
  -d design.mat \
  -t design.con \
  -m ${lesions}/model/lesion_n-map_thr.nii.gz \
  -n 10000 \
  -D \
  -T

# Threshold for significance
...

# Warp from MNI to wmfod space for FBA analysis
mri_synthmorph \
  





