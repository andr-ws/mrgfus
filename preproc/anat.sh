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
raw=${base}/rawdata
der=${base}/derivatives

find ${raw} -type d -name 'sub-*' | sort -V | while read -r dir; do
  sub=$(basename ${dir})
  for ses in ses-01 ses-02 ses-03; do
    pre_anat=${raw}/${sub}/${ses}/anat
    post_anat=${der}/anat/${sub}/${ses}
    mkdir -p ${post_anat}

    # Minimal preprocessing
    for modality in T1w FGATIR; do
      echo "Reorienting and cropping ${sub}..."
      fslreorient2std \
        ${pre_anat}/${sub}_${ses}_acq-${modality}.nii.gz \
        ${post_anat}/${sub}_${ses}_acq-${modality}_reori.nii.gz
	
      robustfov \
        -i ${post_anat}/${sub}_${ses}_acq-${modality}_reori.nii.gz \
        -r ${post_anat}/${sub}_${ses}_acq-${modality}_reori-fov.nii.gz

      echo "Bias-field correcting ${sub}..."
      N4BiasFieldCorrection \
  		-d 3 \
    	-i ${post_anat}/${sub}_${ses}_acq-${modality}_reori-fov.nii.gz \
      -o ${post_anat}/${sub}_${ses}_acq-${modality}_reori-fov-bias.nii.gz

 echo "Brain extracting ${sub}..."
      mri_synthstrip \
        --image ${post_anat}/${sub}_${ses}_acq-${modality}_reori-fov-bias.nii.gz \
        --out ${post_anat}/${sub}_${ses}_acq-${modality}_reori-fov-bias_brain.nii.gz \
        --mask ${post_anat}/${sub}_${ses}_acq-${modality}_reori-fov-bias_brain_mask.nii.gz
    done

    # Co-register FGATIR and T1w modalities
    echo "Co-registering FGATIR to T1w for ${sub}..."
    antsRegistrationSyNQuick.sh \
    	-d 3 \
      -f ${post_anat}/${sub}_${ses}_acq-T1w_brain.nii.gz \
      -m ${post_anat}/${sub}_${ses}_acq-FGATIR_brain.nii.gz \
      -o ${post_anat}/${sub}_${ses}_acq-FGATIR_space-T1w_
    
    echo "Brain extracting ${sub}..."
  done
done
