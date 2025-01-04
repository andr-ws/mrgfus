source=/Volumes/LA_4TB/mrgfus/sourcedata/NEW
raw=/Volumes/LA_4TB/mrgfus/rawdata

for dir in ${source}/sub-*; do
  sub=$(basename ${dir})
  
  mkdir ${source}/${sub}/convert

  /Applications/MRIcroGL.app/Contents/Resources/dcm2niix \
    -o ${raw}/${sub}/ \
    -z y \
    -f ${sub}_%t_%d \
    ${source}/DICOM

done
