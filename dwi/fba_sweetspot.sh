# Code here to perform FBA analysis based on lesion mapping results

mri_synthmorph \
  apply \
  -m nearest /Users/a9ws/imaging/datasets/mrgfus/derivatives/fba/template/study_template/wmfod/MNI-wmfod_template_warp.nii.gz \
  /Users/a9ws/imaging/datasets/mrgfus/derivatives/fba/template/study_template/tracts/nmap_thr.nii.gz \
  /Users/a9ws/imaging/datasets/mrgfus/derivatives/fba/template/study_template/tracts/nmap_thr-wmfod.nii.gz
