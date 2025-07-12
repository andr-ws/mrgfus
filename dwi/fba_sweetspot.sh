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

# Store the metric value of the sweetspot per session (else NA)
for metric in fdc log_fc fd; do
  mean_metric_csv=${fba}/analysis/sweetspot/mean-${metric}.csv
  echo "sub-id,ses-01,ses-02,ses-03" > ${mean_metric_csv}
  mkdir -p ${fba}/analysis/sweetspot/mean_${metric}_maps
done

# Loop through each subject
while read -r sub; do
  echo "Processing $sub"

  for metric in fdc log_fc fd; do
    # Initialize CSV line
    line="$sub"
    
    # Loop through sessions
    for ses in ses-01 ses-02 ses-03; do
      metric_mif=${fba}/template/study_template/metrics/smoothed/${metric}/${sub}_${ses}.mif

      # Check if input file exists
      if [[ -f "$metric_mif" ]]; then
        metric_map=${fba}/analysis/sweetspot/mean_${metric}_maps/${sub}_${ses}.nii.gz

        fixel2voxel \
          ${metric_mif} \
          mean \
          ${metric_map} \
          -weighted ${metric_mif}
          
        fslmaths \
          ${metric_map} \
          -mul \
          ${fba}/analysis/sweetspot/nmap_thr-wmfod.nii.gz \
          ${metric_map}
          
        value=$(fslstats "$metric_map" -M)
        
        newline="${line},${value}"

      else
        echo "  $ses: Missing $metric file."
        newline="${line},NA"
      fi
    done
    echo "$newline" >> ${fba}/analysis/sweetspot/mean-${metric}.csv
  done
done < ${hemis}
