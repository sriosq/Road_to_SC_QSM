# <div align="center">**On the road to QSM of Spinal Cord**</div>

# Theory 

Quantitative Susceptibility Mapping (QSM) is a group of experimental methods that seeks in providing a MR image whose contrast is given by the magnetic suscpetibility of the voxel. Magnetic susceptibility, chi ($\chi$) is the measure of how readily a particle will get magnetized to a external magnetic field. If the particle feels attracted to the field it is said to have a paramagnetic behaviour whereas if it feels repelled to it, it possesses diamagnetic behaviour. There is a third state: Ferromagnetism, but this behaviour is not present in the human body and therefore not quantified for this work.

QSM can be split into image acquisition, image processing and analysis elements [4]. The research focus varies depending of the section of interest. Image acquisition comprises pulse sequences & protocol as well as coil combination, saving and exporting of the data. The processing section is a post-processing pipeline that uses phase information, it begins with Phase Unwrapping and Echo Combination, Mask creation, Background Field Removal and finally Dipole Inversion. 


# Motivation

QSM is a very powerfull tool for non-invasive characterization of iron deposition, estimating voxel oxygenation and geometry, differentiating blood and calcium products as well as studying demnyelinating lesions in white matter. But QSM is still a novel technique, Ithe International Society for Magnetic Resonance in Medicin (ISMRM) has recently published recommendations for implementing QSM for clinical research [4], but it is currently limted to brain analysis. Moreover, most of the parameters used throughout the processing pipeline for QSM have been adjusted for brain images, which underscores the necessity for extending such studies. Here at the NeuroPoly lab in Montreal, one of our primary research focusis is the spinal cord. The motivitation behind this repository is to address the need for extending QSM studies to the spinal cord. In order to perform QSM of the spinal cord, we require a ground truth; therefore the focus is to develop a susceptibility distribution phantom for validating synthetic MRI acquisition. Simulating MRI acquisition will provide means to optimize acquisition parameters for QSM processing pipelines.

# Phantom Creation

This repository uses the publicly available dataset provided by Gatidis et al [2], an annotated Positron Emission Tomography/Computed Tomography (PET/CT) study. The reason CT was chosen instead of a MRI dataset is that with CT data we can use a publicly available tool: Total Segmentator [3]. Total Segmentator will recieve a Nifti Image from CT data and will output a labeled Nifti file where every label corresponds to one of the 117 labels the model was trained to differentiate. This labeled output is of interest as it has the shape of the body and segmentations for different regions such as organs, bones or muscle.

We then group the labels separating them based on the susceptibility value relative to free space that we can use to create a new volume with differentiation of regions based only on the susceptibility. This volume can then be exported as a Nifti file, by default named sus_dist.nii.gz.

We used object oriented python programming to create 2 classes. The first class is the parent class that recieves as input a Nifti file containing the labeled body from Total Segmentator. When an instance of this class is created, using the look up table from Wasserthal et al. [3], all the labels from the file are instanced as Label class (daughter class). When a Label object is created it is automatically attributed relaxation times and proton density values based on literature values depending on the label ID.

An interactive jupyter notebook is provided with comments and an example workflow to understand how the code works and how to use the different implemented methods as the code from this repository is not limited to the usage of the example Nifti file used. It will work with any Nifti file that is an output of Total Segmentator. 

A parcellation color map is provided that can be used with ITK-snap in the file [here](parcellation_itk.txt). This file encodes labels 1 to 48 with a name according to the label name it will have once the code is ran, labels 49 to 67, 72, 73, 77 to 86 are named extra as they don't have a fixed name under the label class.

# Folder Structure

Inside the data folder you will find a raw data folder with the example output of Total Segmentator when using the aforementioned dataset as well as a file containing the relaxation times selected for the assignment on the labels.
 
In the post processing folder you will find a file that is used to recover the body's fat and shape as total segmentator will loose the soft tissue that is arround the organs and that is also limiting the body from the outside. There are also 3 other files: final_total_seg is the file that takes into account the fat and muscle inside the body, new_sc_label file addresses the issue that Total Segmentator will not segment the CSF surrounding the spinal cord and will label as spinal cord the region of the body that is really the Spinal canal and final_sc_seg is the file that introduces 2 new labels to the segmentation: SC_CSF for the cerobrospinal fluid surrounding the spinal cord and Spinal Cord. In order to create this labels, the CT Nifti file was registered to the PAM50 template [5] using Spinal Cord Toolbox [6]. All the necessary files for adding a new label to the segmentation are inside the reslicing folder. This opens an interesting feature: one can add as many labels as desired. When creating and adding a new label it is imperative to select a number greater than 120 so that the label doesn't overlap with the already existing labels from total segmentation. 

The code is inside the simulate_mri folder. The Volume.py file contains the parent class that creates the labels by calling the SegmentationLabel class inside the label.py file. The interactive jupyter notebook will show the usefull methods integrated in the class.

# Data Simulation

The susceptibility distribution phantom (sus_dist.nii.gz) can be used for simulating GRE data acquisition. In the first iteration of the code it is using the following equation:

$$ \text{protonDensity} \cdot \sin(\text{FA}) \cdot \exp\left(-\frac{TE}{T2^*} - \text{sign} \cdot i \cdot \gamma \cdot \delta B_0 \cdot TE\right) $$

The equation is a shorter modified version of the code from the QSM reconstruction challenge 2.0 [7] which assumes sufficient time in between acquisitions (TR) for the longitudinal magnetization to recover to its equilibrium value.

From this equation ($\delta B_0$) is the convolution of the dipole kernel with the susceptibility distribution created earlier.

The code is able to simulate multiple echoes, after the simulation the magnitude and phase can be saved with the appropriate methods.



# References
[1] Chan, K.-S., Marques, J.P., 2021. SEPIA—Susceptibility mapping 
pipeline tool for phase images. Neuroimage 227, 117611. </br>
[2] Gatidis, S., Hepp, T., Früh, M. et al. A whole-body FDG-PET/CT Dataset with manually annotated Tumor Lesions. Sci Data 9, 601 (2022). https://doi.org/10.1038/s41597-022-01718-3 || Attribution to the data set link : https://creativecommons.org/licenses/by/4.0/ </br>
[3] Wasserthal, J., Breit, H.-C., Meyer, M.T., Pradella, M., Hinck, D., Sauter, A.W., Heye, T., Boll, D., Cyriac, J., Yang, S., Bach, M., Segeroth, M., 2023. TotalSegmentator: Robust Segmentation of 104 Anatomic Structures in CT Images. Radiology: Artificial Intelligence. https://doi.org/10.1148/ryai.230024 </br>
[4] QSM Consensus Organization Committee; Bilgic B, Costagli M, Chan KS, Duyn J, Langkammer C, Lee J, Li X, Liu C, Marques JP, Milovic C, Robinson S, Schweser F, Shmueli K, Spincemaille P, Straub S, van Zijl P, Wang Y; ISMRM Electro-Magnetic Tissue Properties Study Group. Recommended Implementation of Quantitative Susceptibility Mapping for Clinical Research in The Brain: A Consensus of the ISMRM Electro-Magnetic Tissue Properties Study Group. ArXiv [Preprint]. 2023 Jul 5:arXiv:2307.02306v1. Update in: Magn Reson Med. 2024 May;91(5):1834-1862. doi: 10.1002/mrm.30006. PMID: 37461418; PMCID: PMC10350101.
[5] Benjamin De Leener, Vladimir S. Fonov, D. Louis Collins, Virginie Callot, Nikola Stikov, Julien Cohen-Adad, PAM50: Unbiased multimodal template of the brainstem and spinal cord aligned with the ICBM152 space, NeuroImage, Volume 165, 2018, Pages 170-179, ISSN 1053-8119, https://doi.org/10.1016/j.neuroimage.2017.10.041. </br>
[6] De Leener B, Levy S, Dupont SM, Fonov VS, Stikov N, Louis Collins D, Callot V, Cohen-Adad J. SCT: Spinal Cord Toolbox, an open-source software for processing spinal cord MRI data. Neuroimage 2017. https://spinalcordtoolbox.com/index.html </br>
[7] Marques JP, Meineke J, Milovic C, et al. QSM reconstruction challenge 2.0: A realistic in silico head phantom for MRI data simulation and evaluation of susceptibility mapping procedures. Magn Reson Med. 2021; 86: 526–542. https://doi.org/10.1002/mrm.28716 </br>