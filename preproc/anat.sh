#!/bin/bash

# Raw T1w/FGATIR handling, minimal pre-processing (reorientation, neck crop, and brain extraction), bias-correction and coreg

export FREESURFER_HOME=/Applications/freesurfer/dev
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

  for ses in ses-preop ses1-postop ses2-postop; do
    
    mkdir -p "${derivatives}/anat/${sub}/${ses}/"

    for modality in T1w FGATIR; do
      echo "Reorienting and cropping ${sub}..."
      fslreorient2std \
        "${rawdata}/${sub}/${ses}/anat/${sub}_acq-${modality}.nii.gz" \
        "${derivatives}/anat/${sub}/${ses}/${sub}_desc-min_proc-${modality}.nii.gz"
	
      robustfov \
        -i "${derivatives}/anat/${sub}/${ses}/${sub}_desc-min_proc-${modality}.nii.gz" \
        -r "${derivatives}/anat/${sub}/${ses}/${sub}_desc-min_proc-${modality}.nii.gz"
    done

    # Co-register FGATIR to T1w
    echo "Co-registering FGATIR to T1w for ${sub}..."
      
    antsRegistrationSyNQuick.sh \
    	-d 3 \
        -f "${derivatives}/anat/${sub}/${ses}/${sub}_desc-min_proc-T1w.nii.gz" \
        -m "${derivatives}/anat/${sub}/${ses}/${sub}_desc-min_proc-FGATIR.nii.gz" \
        -o "${derivatives}/anat/${sub}/${ses}/${sub}_desc-min_proc-FGATIR_space-T1w"

    mv "${derivatives}/anat/${sub}/${ses}/${sub}_desc-min_proc-FGATIR_space-T1wWarped.nii.gz" "${derivatives}/anat/${sub}/${ses}/${sub}_desc-min_proc-FGATIR_space-T1w.nii.gz"
    
    echo "Brain extracting ${sub}..."

    for modality in T1w FGATIR_space-T1w; do
    	mri_synthstrip \
          --image "${derivatives}/anat/${sub}/${ses}/${sub}_desc-min_proc-${modality}.nii.gz" \
          --out "${derivatives}/anat/${sub}/${ses}/${sub}_desc-min_proc-${modality}_brain.nii.gz" \
          --mask "${derivatives}/anat/${sub}/${ses}/${sub}_desc-min_proc-${modality}_brain_mask.nii.gz"
    done
  done
done
