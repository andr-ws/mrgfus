

for dir in ${SUBJECTS_DIR}/sub-*; do
  sub=$(basename ${dir})

  for hemi in lh rh; do
    pathname=${sub}_${ses}.long.sub-001_base
    parc=500.aparc
    for ses in ses-01 ses-02 ses-03; do

      # need to enter code here for transforming fsaverage surface
      # generating annot files etc

      mris_anatomical_stats \
      -a ${SUBJECTS_DIR}/${sub}/long/${pathname}/label/${hemi}.${parc}.annot \
      -f ${SUBJECTS_DIR}/${sub}/long/${pathname}/stats/${hemi}.${parc}.stats  \
      ${sub}/long/${pathname} \
      ${hemi}
    done
  done
done
