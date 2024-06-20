# <div align="center">**On the road to QSM of Spinal Cord**</div>



# Theory 

Quantitative Susceptibility Mapping (QSM) is a group of experimental methods that seeks in providing a MR image whose contrast is given by the magnetic suscpetibility of the voxel. Magnetic susceptibility, chi (\chi) is the measure of how readily a particle will get magnetized to a external magnetic field. If the particle feels attracted to the field it is said to have a paramagnetic behaviour whereas if it feels repelled to it, it possesses diamagnetic behaviour. There is a third state: Ferromagnetism, but this behaviour is not present in the human body and therefore not quantified for this work.

QSM can be split into image acquisition, image processing and analysis elements [4]. The research focus varies depending of the section of interest. Image acquisition comprises pulse sequences & protocol as well as coil combination, saving and exporting of the data. The processing section is a post-processing pipeline that uses phase information, it begins with Phase Unwrapping and Echo Combination, Mask creation, Background Field Removal and finally Dipole Inversion. 


# Motivation

QSM is a very powerfull tool for non-invasive characterization of iron deposition, estimating voxel oxygenation and geometry, differentiating blood and calcium products as well as studying demnyelinating lesions in white matter. But QSM is still a novel technique, Ithe International Society for Magnetic Resonance in Medicin (ISMRM) has recently published recommendations for implementing QSM for clinical research [4], but it is currently limted to brain analysis. Moreover, most of the parameters used throughout the processing pipeline for QSM have been adjusted for brain images, which underscores the necessity for extending such studies. Here at the NeuroPoly lab in Montreal, one of our primary research focusis is the spinal cord. The motivitation behind this repository is to address the need for extending QSM studies to the spinal cord. In order to perform QSM of the spinal cord, we require a ground truth; therefore the focus is to develop a susceptibility distribution phantom for validating synthetic MRI acquisition. Simulating MRI acquisition will provide means to optimize acquisition parameters for QSM processing pipelines.

# Phantom Creation

This repository uses the publicly available dataset provided by Gatidis et al [2], an annotated Positron Emission Tomography/Computed Tomography (PET/CT) study. The reason CT was chosen instead of a MRI dataset is that with CT data we can use a publicly available tool: Total Segmentator [3]. Total Segmentator will recieve a Nifti Image from CT data and will output a labeled Nifti file where every label corresponds to one of the 117 labels the model was trained to differentiate. This labeled output is of interest as it has the shape of the body and segmentations for different regions such as organs, bones or muscle.

We then group the labels separating them based on the susceptibility value relative to free space that we can use to create a new volume with differentiation of regions based only on the susceptibility. This volume can then be exported as a Nifti file, by default named sus_dist.nii.gz.

![alt text](image.png)

We used object oriented python programming to create 2 classes. The first class is the parent class that recieves as input a Nifti file containing the labeled body from Total Segmentator. When an instance of this class is created, given the look up table from Wasserthal et al. [3] we 

####### POLISH THIS LATER ######

Inside the raw data folder you will find the CT whole body Nifti from Gatidis et al. As well as a Nifti file containing the segmentation of the whole body CT done with TotalSegmentator [3].

The post processing folder contains the files used in the interactive notebook. The file from total segmentator contains 114 label. Each label has a susceptibility value. The final segmentation includes a label for fat which is not included in the initial Total Segmentator.

Explain about the way the classes work

Explain about the susceptibility selection criteria

Explain about the susceptibility distribution creation > Acquisition simulation > Importance and applications



For the parcellation color map:
From Label id 1 to 48 everything has a name. 
Then from 49 to 67 everything is automatically labeled "extra"; 72, 73, from 77 to 86 also


Attribution to the data set link : https://creativecommons.org/licenses/by/4.0/ 
It was used to pass through total segmentator CT and then used as template to label and create susceptibility distribution volume. Susceptibility distribution volume is used for simulating MRI acquistion.
# References
[1] Chan, K.-S., Marques, J.P., 2021. SEPIA—Susceptibility mapping 
pipeline tool for phase images. Neuroimage 227, 117611. </br>
[2] Gatidis, S., Hepp, T., Früh, M. et al. A whole-body FDG-PET/CT Dataset with manually annotated Tumor Lesions. Sci Data 9, 601 (2022). https://doi.org/10.1038/s41597-022-01718-3 </br>
[3] Wasserthal, J., Breit, H.-C., Meyer, M.T., Pradella, M., Hinck, D., Sauter, A.W., Heye, T., Boll, D., Cyriac, J., Yang, S., Bach, M., Segeroth, M., 2023. TotalSegmentator: Robust Segmentation of 104 Anatomic Structures in CT Images. Radiology: Artificial Intelligence. https://doi.org/10.1148/ryai.230024 </br>
[4] QSM Consensus Organization Committee; Bilgic B, Costagli M, Chan KS, Duyn J, Langkammer C, Lee J, Li X, Liu C, Marques JP, Milovic C, Robinson S, Schweser F, Shmueli K, Spincemaille P, Straub S, van Zijl P, Wang Y; ISMRM Electro-Magnetic Tissue Properties Study Group. Recommended Implementation of Quantitative Susceptibility Mapping for Clinical Research in The Brain: A Consensus of the ISMRM Electro-Magnetic Tissue Properties Study Group. ArXiv [Preprint]. 2023 Jul 5:arXiv:2307.02306v1. Update in: Magn Reson Med. 2024 May;91(5):1834-1862. doi: 10.1002/mrm.30006. PMID: 37461418; PMCID: PMC10350101.