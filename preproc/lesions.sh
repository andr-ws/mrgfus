#! /bin/bash

# Lesions manually defined
raw=/Volumes/LA_4TB/datasets/mrgfus/rawdata
der=/Volumes/LA_4TB/datasets/mrgfus/derivatives
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

  antsApplyTransforms -d 3 \
    -i ${lesions}/masks/nii/${sub}_lesion.nii.gz \
    -r ${MNI} \
    -t ${der}/anat/${sub}/ses-01/${sub}_T1w-MNI_05mm_1Warp.nii.gz \
    -t ${der}/anat/${sub}/ses-01/${sub}_T1w-MNI_05mm_0GenericAffine.mat \
    -t ${lesions}/data/${sub}/${sub}_immT2w-ses-01_T1w_0GenericAffine.mat \
    -o ${lesions}/masks/mni/${sub}/${sub}_lesion.nii.gz \
    --interpolation NearestNeighbor
done

# Modelling
mkdir ${lesions}/model

# Create an N-mapmrmath \
mrmath \
${lesions}/masks/mni/*/*lesion.nii.gz \
sum \
${lesions}/model/n-map.nii.gz





