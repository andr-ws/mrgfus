#!/bin/bash
#SBATCH --job-name=recon-fs_array
#SBATCH --partition=nodes
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32GB
#SBATCH --time=48:00:00
#SBATCH --array=1-10%5      # 5 is concurrent jobs to process with 10 being the total array size
#SBATCH --mail-type=END
#SBATCH --mail-user=luke.andrews@liverpool.ac.uk

# Code description: runs freesurfer on multi-time-point T1 data
# Organised as ses-preop ses1-postop ses2-postop

# Requirements:
module avail apps/freesurfer/7.4.1
export FREESURFER_HOME=/mnt/data1/users/software/freesurfer/7.4.1
source $FREESURFER_HOME/SetUpFreeSurfer.sh
export SUBJECTS_DIR=/mnt/data1/users/a9ws/freesurfer

# Load subjects from file
subjects_file=/mnt/data1/users/a9ws/subjects_list_fs.txt
mapfile -t subjects < "$subjects_file"

# Determine the subject for this SLURM_ARRAY_TASK_ID
subject_path=${subjects[$SLURM_ARRAY_TASK_ID - 1]}
sub=$(basename "${subject_path}")

# run recon-all
for ses in ses ses-preop ses1-postop ses2-postop; do
  recon-all -s ${sub}_${ses} -i ${SUBJECTS_DIR}/${sub}/T1w-${ses}.nii \
  -all \
  -qcache \
  -openmp 5 # Match to the number of concurrent jobs
done
