#!/bin/bash

export SUBJECTS_DIR=/Volumes/LA_4TB/datasets/mrgfus/derivatives/freesurfer
cd "$SUBJECTS_DIR" || exit 1

for dir in ${SUBJECTS_DIR}/sub-*; do
  sub=$(basename "$dir")
  base="${sub}_base"
  basecmd="recon-all -base $base"
  found=0

  # Create symlinks for each session
  for ses in ses-01 ses-02 ses-03; do
    src="${dir}/${ses}"
    link="${SUBJECTS_DIR}/${sub}_${ses}"
    if [ -d "$src" ] && [ ! -e "$link" ]; then
      ln -s "$src" "$link"
    fi
  done

  # Append valid -tp flags
  for ses in ses-01 ses-02 ses-03; do
    ses_dir="${dir}/${ses}"
    if [ -d "$ses_dir" ]; then
      basecmd+=" -tp ${sub}_${ses}"
      found=1
    fi
  done

  # Run base recon-all
  if [ "$found" -eq 1 ]; then
    basecmd+=" -all"
    echo "Running: $cmd"
    eval "$cmd"
  else
    echo "No valid sessions found for $sub"
  fi

  # run -long on each session
  for ses in ses-01 ses-02 ses-03; do
    if [ -d "${dir}/${ses}" ]; then
      echo "Running longitudinal for ${sub}_${ses} against $base"
      recon-all -long ${sub}_${ses} $base -all
    fi
  done

  # Clean up directories
  run=long.${sub}_base
  mkdir ${dir}/long
  mv ${SUBJECTS_DIR}/${base} ${dir}/long/
  for ses in ses-01 ses-02 ses-03; do
    if [ -d "${dir}/${ses}" ]; then
      mv ${SUBJECTS_DIR}/${sub}.${run} ${dir}/long/
    fi
  done
  
done
