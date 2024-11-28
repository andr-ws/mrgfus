#!/bin/bash

# Raw T1w/FGATIR handling, minimal pre-processing (reorientation, neck crop, and brain extraction), bias-correction and coreg

export FREESURFER_HOME=/Applications/freesurfer
source $FREESURFER_HOME/SetUpFreeSurfer.sh

# Ensure the script is run with at least one argument
if [ "$#" -eq 0 ]; then
    echo "Usage: bash $0 <path/to/dataset>"
    exit 1
fi

# Base path of the dataset
base=$1

# Validate base path
if [ ! -d "$base" ]; then
    echo "Error: The provided base path does not exist."
    exit 1
fi

# Define directories
global="~/imaging/global"
rawdata="${base}/rawdata"
derivatives="${base}/derivatives"

find "${rawdata}" -type d -name 'sub-*' | sort -V | while read -r dir; do
  # Extract subject-id and create directory
  sub=$(basename "${dir}")
  mkdir -p "${derivatives}/anat/${sub}/proc"

  for ses in ses-preop ses1-postop ses2-postop; do

    for modality in T1w FGATIR; do
      echo "Reorienting and cropping ${sub}..."
      fslreorient2std \
        "${rawdata}/${sub}/${ses}/anat/${sub}_acq-${modality}.nii.gz" \
        "${derivatives}/anat/${sub}/proc/${sub}_desc-min_proc_${ses}-${modality}.nii.gz"
	
      robustfov \
        -i "${derivatives}/anat/${sub}/proc/${sub}_desc-min_proc_${ses}-${modality}.nii.gz" \
        -r "${derivatives}/anat/${sub}/proc/${sub}_desc-min_proc_${ses}-${modality}.nii.gz"

      echo "Biasfield correcting ${sub}..."
      N4BiasFieldCorrection \
        -d 3 \
        -i "${derivatives}/anat/${sub}/proc/${sub}_desc-min_proc_${ses}-${modality}.nii.gz" \
        -o "${derivatives}/anat/${sub}/proc/${sub}_desc-bias_cor_${ses}-${modality}.nii.gz"
    done

    for image in min_proc bias_cor; do
      # Co-register T2w to T1w MRI
      echo "Co-registering FGATIR to T1w for ${sub}..."
      antsRegistrationSyNQuick.sh \
        -d 3 \
        -f "${derivatives}/anat/${sub}/proc/${sub}_desc-${image}_${ses}-T1w.nii.gz" \
        -m "${derivatives}/anat/${sub}/proc/${sub}_desc-${image}_${ses}-FGATIR.nii.gz" \
        -o "${derivatives}/anat/${sub}/proc/${sub}_desc-${image}_${ses}-FGATIR_space-T1w"

      mv "${derivatives}/anat/${sub}/proc/${sub}_desc-${image}_${ses}-FGATIR_space-T1wWarped.nii.gz" \
      "${derivatives}/anat/${sub}/proc/${sub}_desc-${image}_${ses}-FGATIR_space-T1w.nii.gz"
    
      echo "Brain extracting ${sub}..."

      for modality in T1w FGATIR; do
        mri_synthstrip \
          --image "${derivatives}/anat/${sub}/proc/${sub}_desc-${image}_${ses}-${modality}.nii.gz" \
          --out "${derivatives}/anat/${sub}/proc/${sub}_desc-${image}_${ses}-${modality}_brain.nii.gz" \
          --mask "${derivatives}/anat/${sub}/proc/${sub}_desc-${image}_${ses}-${modality}_brain_mask.nii.gz"
      done
    done
  done
done
