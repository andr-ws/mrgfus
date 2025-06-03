% MATLAB code to flip rh lesions to lh

ids = readlines('/Volumes/LA_4TB/datasets/mrgfus/derivatives/study_files/rh_subs.txt');

for i = 1:length(ids)
    id = ids{i};
    
    input_file = sprintf('/Users/neuro-239/imaging/datasets/mrgfus/derivatives/lesions/masks/mni/%s/%s_lesion.nii.gz', id, id);
    output_file = sprintf('/Users/neuro-239/imaging/datasets/mrgfus/derivatives/lesions/masks/mni/%s/%s_flesion.nii.gz', id, id);
    
    ea_flip_lr_nonlinear(input_file, output_file, 1);
end
