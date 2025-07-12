#!/bin/bash

# Code here to perform FBA analysis based on lesion mapping results

base=~/imaging/datasets/mrgfus
raw=${base}/rawdata
dwi=${base}/derivatives/data/dwi
fba=${base}/derivatives/fba
sf=${base}/derivatives/study_files

mri_synthmorph \
  apply \
  -m nearest ${fba}/template/study_template/wmfod/MNI-wmfod_template_warp.nii.gz \
  ${fba}/template/study_template/tracts/nmap_thr.nii.gz \
  ${fba}/template/study_template/tracts/nmap_thr-wmfod.nii.gz

tckgen \
  ${fba}/template/study_template/wmfod/wmfod_template.mif \
 ${fba}/template/study_template/tracts/dn-nmap.tck \
  -seed_image ${fba}/template/study_template/tracts/R_DN.mif \
  -select 10000 \
  -include ${fba}/template/study_template/tracts/nmap_thr-wmfod.nii.gz \
  -seed_unidirectional -stop



# Metrics for the nmap GAMM

# Output CSV
combined_csv="${fba}/analysis/sweetspot/combined_metrics_long.csv"
echo "sub-id,timepoint,metric,value" > "${combined_csv}"

# Map session names to timepoint labels
declare -A tp_map
tp_map=( ["ses-01"]="tp1" ["ses-02"]="tp2" ["ses-03"]="tp3" )

# Read subject list
while read -r sub; do
  echo "Processing $sub"

  for metric in fdc log_fc fd; do
    for ses in ses-01 ses-02 ses-03; do

      tp_label=${tp_map[$ses]}
      metric_mif="${fba}/template/study_template/metrics/smoothed/${metric}/${sub}_${ses}.mif"
      metric_map="${fba}/analysis/sweetspot/mean_${metric}_maps/${sub}_${ses}.nii.gz"

      if [[ -f "$metric_mif" ]]; then
        # Compute mean metric in sweet spot
        fixel2voxel "${metric_mif}" mean "${metric_map}" -weighted "${metric_mif}"
        fslmaths "${metric_map}" -mul "${fba}/analysis/sweetspot/nmap_thr-wmfod.nii.gz" "${metric_map}"
        value=$(fslstats "${metric_map}" -M)
        echo "${sub},${tp_label},${metric},${value}" >> "${combined_csv}"
      else
        echo "  $ses: Missing $metric file."
        echo "${sub},${tp_label},${metric},NA" >> "${combined_csv}"
      fi

    done
  done
done < "${hemis}"

