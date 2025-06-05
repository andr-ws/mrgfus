#! /bin/bash

# Lesions manually defined
raw=/Volumes/LA_4TB/datasets/mrgfus/rawdata
der=/Volumes/LA_4TB/datasets/mrgfus/derivatives
sf=${der}/study_files
lesions=${der}/lesions
MNI=/Applications/leaddbs/templates/space/MNI152NLin2009bAsym/t1_brain.nii.gz

# Need to generate coreg btwn lesion T2 and ses-01 T1
for dir in ${raw}/sub-*; do
  sub=$(basename ${dir})
  mkdir -p ${lesions}/data/${sub}
  T2=${dir}/imm_postop/${sub}_acq-T2w.nii.gz
  T1=${der}/anat/${sub}/ses-01/${sub}_ses-01_acq-T1w_biasco.nii.gz

  # Rigid affine transform
  antsRegistrationSyNQuick.sh \
    -d 3 \
    -f ${T1} \
    -m ${T2} \
    -o ${lesions}/data/${sub}/${sub}_immT2w-ses-01_T1w_ \
    -t r \
    -n 12
  
  # Warp the lesion masks
  mkdir -p ${lesions}/masks/mni/${sub}

  antsApplyTransforms \
    -d 3 \
    -i ${lesions}/masks/nii/${sub}_lesion.nii.gz \
    -r ${MNI} \
    -t ${der}/anat/${sub}/ses-01/${sub}_T1w-MNI_05mm_1Warp.nii.gz \
    -t ${der}/anat/${sub}/ses-01/${sub}_T1w-MNI_05mm_0GenericAffine.mat \
    -t ${lesions}/data/${sub}/${sub}_immT2w-ses-01_T1w_0GenericAffine.mat \
    -o ${lesions}/masks/mni/${sub}/${sub}_lesion.nii.gz \
    --interpolation NearestNeighbor
done

# Flip rh lesions to lh
matlab \
  -nodisplay \
  -nosplash \
  -r "run('/Users/neuro-239/scripts/flip_rh_lesions.m'); exit;"

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
