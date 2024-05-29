% Go to the folder and select the path for the nii gz :)

path_to_nii = 'C:/Users/User/QSM_data/patient1/fmap_thr_50.nii.gz';
path_to_local_field = 'Sepia_localfield.nii.gz';
path_to_chi_map = 'Sepia_Chimap.nii.gz';
path_fieldmap = 'Sepia_fieldmap.nii.gz';

nii_img = load_untouch_nii(path_to_chi_map);
img_data = nii_img.img;

% Get the number of slices
num_slices = size(img_data, 3);

%%

% Prompt the user for a slice number
slice_number = input('Enter the slice number to display: ');

% Check if the input slice number is valid
if slice_number < 1 || slice_number > num_slices
    error('Invalid slice number. Please enter a number between 1 and %d.', num_slices);
end

% Display the selected slice
figure;
imagesc(img_data(:,:,slice_number));
colormap(gray); 
axis off; % Turn off axis labels
title(['Slice ', num2str(slice_number)]);

% Adjust limits after seeing the complete image
% For spinal cord square box usually is 170,220 (Try iterating though)
x_min = 170;
x_max = 220;
y_min = 155;
y_max = 220;                    
width = x_max - x_min;
height = y_max - y_min;
rectangle('Position', [x_min, y_min, width, height], 'EdgeColor', 'r', 'LineWidth', 2);

xlim([100, 250]);
ylim([150, 250]);

% Adjust contrast using clim
clim([-0.00011283 0.00012595]);
set(gca, 'CLim', clim);

colormap;
colorbar;
h = colorbar;
title(h, 'Hz');
