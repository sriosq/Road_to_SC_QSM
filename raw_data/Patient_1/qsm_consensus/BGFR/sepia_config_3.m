% This file is generated by SEPIA version: v1.2.2.6
% add general Path
sepia_addpath;

% Input/Output filenames
input = struct();
input(1).name = 'H:\My Drive\Maitrise\Project\Data\first qsm trial\Patient_1\qsm_consensus\Phase_unwrapp\Sepia_fieldmap.nii.gz' ;
input(2).name = '' ;
input(3).name = 'H:\My Drive\Maitrise\Project\Data\first qsm trial\Patient_1\qsm_consensus\Phase_unwrapp\Sepia_noisesd.nii.gz' ;
input(4).name = 'H:\My Drive\Maitrise\Project\Data\first qsm trial\sepia_header.mat' ;
output_basename = 'H:\My Drive\Maitrise\Project\Data\first qsm trial\Patient_1\qsm_consensus\BGFR\Sepia' ;
mask_filename = ['H:\My Drive\Maitrise\Project\Data\first qsm trial\Patient_1\qsm_consensus\Phase_unwrapp\Sepia_mask_localfield.nii.gz'] ;

% General algorithm parameters
algorParam = struct();
algorParam.general.isBET = 0 ;
algorParam.general.isInvert = 0 ;
algorParam.general.isRefineBrainMask = 0 ;
% Background field removal algorithm parameters
algorParam.bfr.refine_method = 'None' ;
algorParam.bfr.refine_order = 4 ;
algorParam.bfr.erode_radius = 0 ;
algorParam.bfr.erode_before_radius = 1 ;
algorParam.bfr.method = 'VSHARP' ;
algorParam.bfr.radius = [8:-1:1] ;

sepiaIO(input,output_basename,mask_filename,algorParam);