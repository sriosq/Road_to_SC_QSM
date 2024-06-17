# <div align="center">**Simulating MRI**</div>

# Susceptibility Phantom of the Spinal Cord
Repository for the Brainhack project - summer 2024

# Background

MRI provides a range of image contrasts typically qualitative i.e. signals reflect 
tissue changes but are not direct measures of the tissue properties.

This repository can create a digital phantom of the spinal cord that allows for 
simulation of MR signal acquisition to assess the parameters for calculating 
susceptibility values of the spinal cord.

This repository uses the publicly available dataset provided by Gatidis et al [2].

####### POLISH THIS LATER ######

Inside the raw data folder you will find the CT whole body Nifti from Gatidis et al. As well as a Nifti file containing the segmentation of the whole body CT done with TotalSegmentator [3].

The post processing folder contains the files used in the interactive notebook. The file from total segmentator contains 114 label. Each label has a susceptibility value. The final segmentation includes a label for fat which is not included in the initial Total Segmentator.

Explain about the way the classes work

Explain about the susceptibility selection criteria

Explain about the susceptibility distribution creation > Acquisition simulation > Importance and applications


Attribution to the data set link : https://creativecommons.org/licenses/by/4.0/ 
It was used to pass through total segmentator CT and then used as template to label and create susceptibility distribution volume. Susceptibility distribution volume is used for simulating MRI acquistion.
# References
[1] Chan, K.-S., Marques, J.P., 2021. SEPIA—Susceptibility mapping 
pipeline tool for phase images. Neuroimage 227, 117611. </br>
[2] Gatidis, S., Hepp, T., Früh, M. et al. A whole-body FDG-PET/CT Dataset with manually annotated Tumor Lesions. Sci Data 9, 601 (2022). https://doi.org/10.1038/s41597-022-01718-3 </br>
[3] Wasserthal, J., Breit, H.-C., Meyer, M.T., Pradella, M., Hinck, D., Sauter, A.W., Heye, T., Boll, D., Cyriac, J., Yang, S., Bach, M., Segeroth, M., 2023. TotalSegmentator: Robust Segmentation of 104 Anatomic Structures in CT Images. Radiology: Artificial Intelligence. https://doi.org/10.1148/ryai.230024 </br>
