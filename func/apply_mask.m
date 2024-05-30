% Function to apply a mask



% Load the brain image and its header
% Choosing path:
path_img = 'Sepia_Chimap.nii.gz';
path_mask = 'diluted_sc_mask.nii';


%% 
brain_nii = load_untouch_nii(path_img);
brain_image = double(brain_nii.img);
brain_image(isnan(brain_image)) = 0;
% Load the mask image and its header
mask_nii = load_untouch_nii(path_mask);
mask_image = double(mask_nii.img);

% Apply the mask to the brain image
masked_image = brain_image .* mask_image;
masked_image(isnan(masked_image)) = 0;
% Create a new NIfTI structure for the masked image with the same header as the brain image
masked_nii = mask_nii;
masked_nii.img = masked_image;

% Save the masked image as a new NIfTI file
save_untouch_nii(masked_nii, 'image_masked.nii');
%%
% Load the brain image and its header
brain_nii = load_untouch_nii('Sepia_Chimap.nii.gz');
brain_image = double(brain_nii.img); % Convert to double
% Load the mask image and its header
mask_nii = load_untouch_nii('diluted_sc_mask.nii');
mask_image = double(mask_nii.img); % Convert to double
% Apply the mask to the brain image
masked_image = brain_image .* mask_image;
% Create a new NIfTI structure for the masked image with the same header as the brain image
masked_nii = brain_nii; % Use the brain image header for consistency
masked_nii.img = masked_image;
% Save the masked image as a new NIfTI file
save_untouch_nii(masked_nii, 'Sepia_Chimap_masked.nii')