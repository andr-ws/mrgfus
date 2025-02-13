#!/bin/bash

# Code to compute scanner-specific white matter FODs and tractograms using the synthetic b0 image for ACT and parcellation.

der=~/imaging/datasets/mrgfus/derivatives
fba=${der}/projects/fba
sf=${der}/projects/study_files
  
for dir in ${fba}/data/sub-*; do
  sub=$(basename ${dir})

  for ses in ses-01 ses-02 ses-03; do
    # Setup ACT
    fba_sdir=${dir}/${ses}
    
    mkdir ${fba_sdir}/ACT
    
    # 5tt-generation
    5ttgen fsl ${der}/data/anat/${sub}/${ses}/${sub}_${ses}_acq-T1w_brain.nii.gz ${fba_sdir}/ACT/${sub}_5tt.mif -nocrop -premasked -force

    # Generate interface
    5tt2gmwmi ${fba_sdir}/ACT/${sub}_5tt.mif ${fba_sdir}/ACT/${sub}_gmwmi.mif -force

    # Compute full tractogram
    tckgen -angle 22.5 -maxlen 250 -minlen 10 -power 1.0 -cutoff 0.10 -select 2000000 \
    ${fba_sdir}/fod/${sub}_wmfod.mif \
    -seed_gmwmi ${fba_sdir}/ACT/${sub}_gmwmi.mif \
    -act ${fba_sdir}/ACT/${sub}_5tt.mif \
    ${fba_sdir}/ACT/${sub}.tck \
    -force

    # SIFT and extract weights
    tcksift2 ${fba_sdir}/ACT/${sub}.tck ${fba_sdir}/fod/${sub}_wmfod.mif ${fba_sdir}/ACT/${sub}_weights.txt \
    -act ${fba_sdir}/ACT/${sub}_5tt.mif -fd_scale_gm -out_mu ${fba_sdir}/ACT/${sub}_pc_mu.txt -force

    # Reformat schaefer parcellation
    labelconvert \
    ${der}/data/freesurfer/${sub}/${ses}/parcellations/${sub}_Schaefer2018_400Parcels_7Networks_space-dwi.nii.gz \
    ${sf}/parcellations/lut/mrtrix/Schaefer2018_400Parcels_7Networks_order_mrtrix.txt \
    ${sf}/parcellations/lut/mrtrix/Schaefer2018_400Parcels_7Networks_order_mrtrix_reorder.txt \
    ${fba_sdir}/ACT/Schaefer-400_parc.mif

    # Compute connectome
    tck2connectome ${fba_sdir}/ACT/${sub}.tck ${fba_sdir}/ACT/Schaefer-400_parc.mif ${fba_sdir}/ACT/Schaefer_tck_400.csv \
    -tck_weights_in ${fba_sdir}/ACT/${sub}_weights.txt \
    -symmetric -zero_diagonal
      
    # Delete tractogram
    rm ${fba_sdir}/ACT/${sub}.tck
  done
done
