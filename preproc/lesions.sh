#! /bin/bash

# Lesions manually defined
raw=/Volumes/LA_4TB/datasets/mrgfus/rawdata
der=/Volumes/LA_4TB/datasets/mrgfus/derivatives
lesions=${der}/lesions

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
    
done
