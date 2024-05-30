
clear

path = 'D:/Poly_MSc_Code/libraries_and_toolboxes/sepia';
path2 = 'D:/Poly_MSc_Code/libraries_and_toolboxes/toolboxes';

data_path = 'H:/My Drive/Maitrise/Project/Data/first qsm trial';
addpath(path);
addpath(genpath(path2));
cd(path);

%Run to open sepia 
sepia
%% 
path_to_nii = 'C:/Users/User/QSM_data/patient1/fmap_thr_50.nii.gz';
path_to_local_field = 'Sepia_localfield.nii.gz';
path_to_chi_map = 'Sepia_Chimap.nii.gz';
path_fieldmap = 'Sepia_fieldmap.nii.gz';

nii_img = load_untouch_nii(path_to_chi_map);

img_data = nii_img.img;

% Get the number of slices
num_slices = size(img_data, 3);

% Assuming every slice has the same x y dimension limitss
[rows, cols] = size(img_data(:,:,1));
figure;
% Display each slice separately
for slice = 1:num_slices
    subplot(4, 4, slice); % Change subplot dimensions as needed
    imagesc(img_data(:,:,slice));
    colormap(gray); 
    axis off; % Turn off axis labels
    title(['Slice ', num2str(slice)]);
    % Adjust limits after seeing the complete image
    % For spinal cord square box usually is 170,220 (Try iterating tho)
    x_min = 170;
    x_max = 220;
    y_min = 155;
    y_max = 220;                    
    width = x_max - x_min;
    height = y_max - y_min;
    rectangle('Position', [x_min, y_min, width, height], 'EdgeColor', 'r', 'LineWidth', 2);

    xlim([100, 250]);
    ylim([150, 250]);
    % With clim we adjust the contrast
    % Using contrast inspector from ITK snap we can quickly get this range
    % of values
    clim([-0.00011283 0.00012595]);

    colormap;
    colorbar;
    h = colorbar;
    % If path is path_to_nii4
    title(h,'Hz')
    % else
    %title(h,'Hz');
end