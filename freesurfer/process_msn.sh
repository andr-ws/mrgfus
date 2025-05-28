#! /bin/bash

# Code to resample fsaverage 500 ROI DKT atlas onto longitudinal pipeline surfaces

for dir in ${SUBJECTS_DIR}/sub-*; do
  sub=$(basename ${dir})

  for hemi in lh rh; do
    parc=500.aparc
    
    for ses in ses-01 ses-02 ses-03; do
      pathname=${sub}_${ses}.long.${sub}_base
      
      mri_surf2surf \
        --hemi ${hemi} \
        --srcsubject fsaverage \
        --trgsubject ${sub}/long/${pathname} \
        --sval-annot ${SUBJECTS_DIR}/fsaverage/label/${hemi}.${parc}.annot \
        --tval ${sub}/long/${pathname}/label/${hemi}.${parc}.annot
    
      mris_anatomical_stats \
        -a ${SUBJECTS_DIR}/${sub}/long/${pathname}/label/${hemi}.${parc}.annot \
        -f ${SUBJECTS_DIR}/${sub}/long/${pathname}/stats/${hemi}.${parc}.stats  \
        ${sub}/long/${pathname} \
        ${hemi}
    done
  done
done
