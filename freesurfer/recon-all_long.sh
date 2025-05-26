


export SUBJECTS_DIR=/Volumes/LA_4TB/datasets/mrgfus/derivatives/freesurfer

# Pipeline expects sessions to dirs nested within SUBJECTS_DIR so symlink

cd /Volumes/LA_4TB/datasets/mrgfus/derivatives/freesurfer

ln -s ${SUBJECTS_DIR}/sub-001/ses-01 ${SUBJECTS_DIR}/sub-001_ses-01
ln -s ${SUBJECTS_DIR}/sub-001/ses-02 ${SUBJECTS_DIR}/sub-001_ses-02
ln -s ${SUBJECTS_DIR}/sub-001/ses-03 ${SUBJECTS_DIR}/sub-001_ses-03

recon-all -base sub-001_6m_base \
  -tp sub-001_ses-01 \
  -tp sub-001_ses-02 \
  -all

recon-all -long sub-001_ses-01 sub-001_6m_base -all
recon-all -long sub-001_ses-02 sub-001_6m_base -all

# Repeat for 12m
recon-all -base sub-001_12m_base \
  -tp sub-001_ses-01 \
  -tp sub-001_ses-03 \
  -all

