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

  for ses in ses-01 ses-02 ses-03; do
    
    mkdir -p "${derivatives}/anat/${sub}/${ses}/"

    for modality in T1w FGATIR; do
      
      echo "Reorienting and cropping ${sub}..."
      fslreorient2std \
        "${rawdata}/${sub}/${ses}/anat/${sub}_${ses}_acq-${modality}.nii.gz" \
        "${derivatives}/anat/${sub}/${ses}/${sub}_${ses}_acq-${modality}_reoriented.nii.gz"
	
      robustfov \
        -i "${derivatives}/anat/${sub}/${ses}/${sub}_${ses}_acq-${modality}_reoriented.nii.gz" \
        -r "${derivatives}/anat/${sub}/${ses}/${sub}_${ses}_acq-${modality}_fov.nii.gz"

      echo "Bias-field correcting ${sub}..."

      N4BiasFieldCorrection \
  		-d 3 \
    	-i "${derivatives}/anat/${sub}/${ses}/${sub}_${ses}_acq-${modality}_fov.nii.gz" \
      -o "${derivatives}/anat/${sub}/${ses}/${sub}_${ses}_acq-${modality}_biasco.nii.gz"

 	    echo "Brain extracting ${sub}..."
  	   
      mri_synthstrip \
        --image "${derivatives}/anat/${sub}/${ses}/${sub}_${ses}_acq-${modality}_biasco.nii.gz" \
        --out "${derivatives}/anat/${sub}/${ses}/${sub}_${ses}_acq-${modality}_brain.nii.gz" \
        --mask "${derivatives}/anat/${sub}/${ses}/${sub}_${ses}_acq-${modality}_brain_mask.nii.gz"
 
    done

    # Co-register FGATIR to T1w
    echo "Co-registering FGATIR to T1w for ${sub}..."
      
    antsRegistrationSyNQuick.sh \
    	-d 3 \
      -f "${derivatives}/anat/${sub}/${ses}/${sub}_${ses}_acq-T1w_brain.nii.gz" \
      -m "${derivatives}/anat/${sub}/${ses}/${sub}_${ses}_acq-FGATIR_brain.nii.gz" \
      -o "${derivatives}/anat/${sub}/${ses}/${sub}_${ses}_acq-FGATIR_space-T1w_"

    mv "${derivatives}/anat/${sub}/${ses}/${sub}_${ses}_acq-FGATIR_space-T1w_Warped.nii.gz" "${derivatives}/anat/${sub}/${ses}/${sub}_${ses}_acq-FGATIR_coreg.nii.gz"
    
    echo "Brain extracting ${sub}..."

    rm -r ${derivatives}/anat/${sub}/${ses}/*warp*
    rm -r ${derivatives}/anat/${sub}/${ses}/*.mat*
  
  done
done
