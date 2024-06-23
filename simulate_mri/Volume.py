#Dependencies
import numpy as np
from simulate_mri.label import SegmentationLabel
import nibabel as nib
from simulate_mri import *
from simulate_mri.utils.simulation_functions import create_dipole_kernel, generate_signal, show_slices, optimized_signal
import scipy.ndimage
from simulate_mri.utils.get_dic_values import to_csv_sus, to_csv_relax
import os

# Parent class for the creation of a non-finite biomechanical model of the body
class Volume:
    
    def __init__(self, volume):
        # In this version we correct that the output should be the nifti image
        # This way we can attribute the information from nifti files to the class

        self.nifti = volume # This is now directing to a Nifti file
        self.volume = self.nifti.get_fdata()
        self.dimensions = np.array(self.volume.shape) # It is initially a tuple, but it needs to be an array
        self.uniq_labels = np.unique(self.volume)
        self.segmentation_labels = {} 
        self.sus_dist = np.zeros(self.dimensions)
        self.t2star_vol = np.zeros(self.dimensions)
        self.pd_dist = np.zeros(self.dimensions)
        self.deltaB0 = np.zeros(self.dimensions)
        # The dictionary has keys for every id number and each value 
        # is the corresponding SegmentationLabel daughter class

        self.create_labels()
        # This is the Convention Dictionary for labels_id - names - sus_values

        # Creating folders for the code
        if not os.path.exists("output"):
            os.makedirs('output')
        if not os.path.exists('simulation'):
            os.makedirs('simulation')

        self.magnitude = None
        self.phase = None
        self.real = None
        self.imaginary = None

    def create_labels(self):
       for label_id in self.uniq_labels:
           self.segmentation_labels[label_id] = SegmentationLabel(label_id)

    def create_segmentation_labels(self):

        # The most important labels are:
        # lungs Bone Soft Tissue SpinalCord CSF
        # These are the regions of the body that impact the most
        # For susceptibility values: https://pfeifer.phas.ubc.ca/refbase/files/Truong-MRI-2002-20-759.pdf
        # from bone, fat and softtissue
        # We can create an already pre-defined from convention of label_ids

        # In TotalSegmentator this is labeled Spinal Cord, but
        # it really is the Spinal Canal = Spinal Cord + CSF
        # For SC is (GM + WM)/2 = -9.055
        # From Eva's code GM = -9.03 and WM = -9.08
        self.set_label_name(76,"SpinalCanal")
        self.set_label_susceptibility(76,-9.055)
        # IN this branch of the code, we have added labels for the SC
        # We correct the fact that in Total Seg it was Spinal Canal
        self.set_label_name(256, "spinal_cord")
        self.set_label_susceptibility(256,-9.055)
        self.set_label_name(289, "sc_csf")
        self.set_label_susceptibility(289,-9.05)


        # For the lungs we have [9,10,11,12,13]
        for i in [9,10,11,12,13]:
            self.set_label_name(i,'lungs')
            self.set_label_susceptibility(i,0.4)

        # For the bones we have a lot more labels
        # Vertebrae and ribs. But all of them can have the same value: -11.1
        # Vertebrae list goes from Sacrum to C1
        vertebra_list = np.arange(22,48)
        rib_list = np.append(np.arange(89,114),[68,69,70,71,74,75,88])

        # It was -11.1, but now it is -9.0
        for i in vertebra_list:
            self.set_label_name(i,"bone")
            self.set_label_susceptibility(i,-9.0)

        for i in rib_list:
            self.set_label_name(i, "bone")
            self.set_label_susceptibility(i, -9.0)

        # Last but not least is to give susceptibility values to the organs
        # and soft tissue => susceptibility value of water = -9.05

        sus_water_list = [1, 2, 3, 4, 5, 6, 7, 8, 14, 16]

        self.set_label_name(1, "spleen")
        self.set_label_name(2, "kidney") #Right
        self.set_label_name(3, "kidney") # Left
        self.set_label_name(4, "organ") # Gallblader
        self.set_label_name(5, "liver")
        self.set_label_name(6, "organ") #Stomach
        self.set_label_name(7, "gland") # AdrenalGland
        self.set_label_name(8, "gland") # AdrenalGland
        self.set_label_name(14, "esophagus")
        self.set_label_name(15, "trachea")
        self.set_label_name(16, "gland") # Thyroid


        # Susceptibility value of fat = -8.39

        for i in sus_water_list:
            self.set_label_susceptibility(i, -9.05)

        # For trachea, it should have a susceptibilty value closer to air
        self.set_label_susceptibility(15, 0.4)

        # Soft tissue == water
        # Inside of body == fat


        self.set_label_name(0,"air") # Outside of brain
        self.set_label_susceptibility(0,0.35)

        # If label has not been set it can be considered as fat
        # susceptibility of fat label is considering fat and muscle proportion of the body
        # sus of fat is -7.5, muscle is -9.05 and soft tissue is -9.5
        # Assuming that body is 20% fat, 40% muscle and 40% soft tissue
        # The weighted average of the fat label should be: -8.92

        self.set_label_name(264,"fat")
        self.set_label_susceptibility(264,-8.92)

        # For simulating GRE acquisition we need to set name of organs
        # brain, liver, spleen, kidney
        self.set_label_name(87, "brain")
        self.set_label_name(48,"heart")


        # The other labels missing without names follow

        intestines = [17,18,19,20,21,22,23]
        # small intestine duodenum colon urinary bladder prostate kidney cyst left and right
        for i in intestines:
            self.set_label_name(i,"organ")
            self.set_label_susceptibility(i,-9.05)

        # Lastly if everything was labeled properly what's left are veins and muscles
        # The susceptibility of csf and water is -9.05, of muscle is -9.03 we can combine them
        for i in self.uniq_labels:
            if self.segmentation_labels[i].name == None:
                self.set_label_name(i,"extra")
                self.set_label_susceptibility(i,-9.04)

        # This way we should have everything labeled
    def check_labels(self):
        for i in self.uniq_labels:
            if self.segmentation_labels[i].name == None:
                print("Label: ",self.segmentation_labels[i]["name"]," doesn't have name assigned")

    def set_label_name(self, label_id, name):
        if label_id in self.uniq_labels:
            self.segmentation_labels[label_id].set_name(name)
        else: print(f"Label ID {label_id} not found.")

    def set_label_susceptibility(self, label_id, susceptibility):
        if label_id in self.uniq_labels:
            # SImilar to set_label_name
            self.segmentation_labels[label_id].set_susceptibility(susceptibility)
        else: print(f"Label ID {label_id} not found.")
    def set_label_pd(self,label_id,pd):
        if label_id in self.uniq_labels:
            self.segmentation_labels[label_id].set_pd_val(pd)
        else: print(f"Label ID {label_id} not found.")

    def set_T2star(self, label_id, t2star):
        if label_id in self.uniq_labels:
        # SImilar to set_label_name
            self.segmentation_labels[label_id].set_t2star_val(t2star)
        else:
            print(f"Label ID {label_id} not found.")
 
    def display(self):
        show_slices(self.volume)


    def manual_label(self,id,name,sus):
        if id in self.uniq_labels:
            label = self.segmentation_labels[id]
            label.name = name
            label.sus = sus

    def show_labels(self):
        for i in self.segmentation_labels:
            label = self.segmentation_labels[i]
            print(label) # Calling __str__ from label


    def create_sus_dist(self):
        # Code for create a susceptibility distribution volume
        # Using the label class
        for i in range(self.dimensions[0]):
            for j in range(self.dimensions[1]):
                for k in range(self.dimensions[2]):

                    pixel = self.volume[i,j,k]
                    label = self.segmentation_labels[pixel]
                    suscep = label.susceptibility

                    if suscep == None:
                        # THis means the label is not defined
                        # The only not defined labels are organs
                        # We can consider the suscp of water
                        self.sus_dist[i,j,k] = -9.05
                    else:
                        self.sus_dist[i,j,k] = suscep

        return self.sus_dist

    def save_sus_dist_nii(self):
        # Method to save the susceptibility distribution created to nifti
        temp_img = nib.Nifti1Image(self.sus_dist, affine=self.nifti.affine)
        path = os.path.join('output','sus_dist.nii.gz')
        # Save the new NIfTI image to a file
        nib.save(temp_img,path)
        del temp_img
        del path

    def create_pd_vol(self):
        # This method will use the lookup table of PD values to create a new volume
        # This new volume will use the labels to quickly create a volume with ProtonDensity values

        for i in range(self.dimensions[0]):
            for j in range(self.dimensions[1]):
                for k in range(self.dimensions[2]):

                    pixel = self.volume[i,j,k]
                    label = self.segmentation_labels[pixel]
                    pd = label.PD_val
                    if pd == None:
                        # THis means the label does not have PD defined
                        self.pd_dist[i,j,k] = 0
                    else:
                        # If the label has PD value it will put this value on the volume
                        self.pd_dist[i,j,k] = pd


    def create_t2_star_vol(self):
        # This method will use the lookup table of T2 star values to create a new volume
        # This new volume will use the labels to quickly create a volume with relaxation time

        for i in range(self.dimensions[0]):
            for j in range(self.dimensions[1]):
                for k in range(self.dimensions[2]):

                    pixel = self.volume[i,j,k]
                    label = self.segmentation_labels[pixel]
                    t2star = label.T2star_val
                    if t2star == None:
                        # THis means the label does not have T2 star value defined
                        self.pd_dist[i,j,k] = 0.001
                    else:
                        # If the label has value it will put this value on the volume
                        self.t2star_vol[i,j,k] = t2star


    def save_pd_dist(self):
        # Method to save the proton density distribution created to nifti
        temp_img = nib.Nifti1Image(self.pd_dist, affine=self.nifti.affine)
        # Save the new NIfTI image to a file
        path = os.path.join('output', 'pd_dist.nii.gz')
        nib.save(temp_img,path)
        del temp_img
        del path

    def calculate_deltaB0(self,B0_dir =[0,0,1], T=3):

        voxel_size = self.nifti.header["pixdim"][1:4]
        buffer = 1
        D = create_dipole_kernel(B0_dir,voxel_size,self.dimensions,buffer=buffer)
        padded_dims = self.dimensions*buffer

        fft_chi = np.fft.fftn(self.sus_dist,padded_dims)

        Bz_fft = fft_chi*D

        vol_no_buff = np.real(np.fft.ifftn(Bz_fft))
        # This Bz_fft will come in ppm
        self.deltaB0 = vol_no_buff[0:self.dimensions[0], 0:self.dimensions[1], 0:self.dimensions[2]]
        # To set a Tesla value it will be an input, but for now is default at 3T
        gamma = 42.58e6 # At 1 tesla, will be scaled by T, default 3
        self.deltaB0 = self.deltaB0*T*gamma*1e-6

    def save_deltaB0(self):
        temp_img = nib.Nifti1Image(self.deltaB0, affine=self.nifti.affine)
        path = os.path.join('output', 'freq_map.nii.gz')
        nib.save(temp_img,path)
        del temp_img
        del path
    # This version of the code assumes that TR is long enough for all Longitudinal Magnetization to return
    # to its equilibrium value
    def simulate_measurement(self,FA,TE,B0=3):
        #FA : flip angle
        #T2 star in seconds
        #B0 in Tesla
        # Gamma in rad*Hz/Tesla
        # handedness => % Siemens & Canon = 'left', GE & Philips = 'right'
        # TE should be a list, so we create a new volume
        num_TE = len(TE)
        newVol_dims = list(self.dimensions)
        newVol_dims.append(num_TE)
        # This way we can iterate over the last dimension (TEs)

        self.magnitude = np.zeros(newVol_dims)
        self.phase = np.zeros(newVol_dims)

        gamma = 267.52218744e6*3/B0# Using gamma for 3 Tesla, B0 is optional to change => This is rad*Hz/Tesla
        handedness = 'left'

        for te in range(num_TE):
            for i in range(self.dimensions[0]):
                for j in range(self.dimensions[1]):
                    for k in range(self.dimensions[2]):

                        pixel = self.volume[i,j,k]
                        deltaB0 = self.deltaB0[i,j,k]
                        label = self.segmentation_labels[pixel]
                        pd = label.PD_val
                        t2star = label.T2star_val
                        mag,phase = generate_signal(pd,t2star,FA,te,deltaB0,gamma,handedness)
                        self.magnitude[i,j,k,te] = mag
                        self.phase[i,j,k,te] = phase

    # Line to implement from MATLAB
    # signal = pd*sind(FA)*exp(-TE./T2star-sign*1i*gamma*deltaB0*TE)

    def optimize_measurement(self,FA,TE,B0=3):
        # This code seeks to accomplish the same as the above method but
        # We are trying to optimize by using volumes
        self.create_pd_vol()
        self.create_t2_star_vol()
        self.calculate_deltaB0()

        # TE should be a list, so we create a new volume
        num_TE = len(TE)
        newVol_dims = list(self.dimensions)
        newVol_dims.append(num_TE)
        # This way we can iterate over the last dimension (TEs)

        self.magnitude = np.zeros(newVol_dims)
        self.phase = np.zeros(newVol_dims)

        gamma = 267.52218744e6 * 3 / B0  # Using gamma for 3 Tesla, B0 is optional to change => This is rad*Hz/Tesla
        handedness = 'left'

        for te_idx, TE in enumerate(TE):
            mag,phase = optimized_signal(self.pd_dist, self.t2star_vol, FA, TE, self.deltaB0, gamma, handedness)
            self.magnitude[...,te_idx] = mag
            self.phase[...,te_idx] = phase


    def get_Magnitude(self):
        temp_img = nib.Nifti1Image(self.magnitude, affine=self.nifti.affine)
        path = os.path.join('simulation','magnitude.nii.gz')
        # Save the new NIfTI image to a file
        nib.save(temp_img, path)
        del temp_img
        del path

    def get_Phase(self):
        temp_img = nib.Nifti1Image(self.phase, affine=self.nifti.affine)
        path = os.path.join('simulation','phase.nii.gz')
        # Save the new NIfTI image to a file
        nib.save(temp_img, path)
        del temp_img
        del path
    def get_Real(self):
        temp_img = nib.Nifti1Image(self.real, affine=self.nifti.affine)
        path = os.path.join('simulation','real.nii.gz')
        # Save the new NIfTI image to a file
        nib.save(temp_img, path)
        del temp_img
        del path

    def get_Imaginary(self):
        temp_img = nib.Nifti1Image(self.imag, affine=self.nifti.affine)
        path = os.join('simulation','imaginary.nii.gz')
        # Save the new NIfTI image to a file
        nib.save(temp_img, path)
        del temp_img



    # Implementation of the code to simulate MRI data acquisition
    # This code is more complex as it assumes that there is still Longitudinal Magnetizatation - T1 T2 and M0
    # Values are necessary
    def simulate_signal_hard(self,TE,TR,theta,B0):
        # Theta is a fixed angle // TE can be multiple echo time array or 1 echo time
        # TR is the repetition time

        B0_dir = [0, 0, 1]
        voxel_size = self.nift.header['pixdim'][1:4] # Resolution - voxel size

        # D = create_dipole_kernel(B0_dir,voxel_size,self.dimensions) # Now its an attribute

        # This variable is temporable and will later be deleted
        chitemp = np.ones([2*d for d in self.dimensions]) *self.sus_dist[-1,-1,-1]

        chitemp[:self.dimensions[0], :self.dimensions[1], :self.dimensions[2]] = self.sus_dist

        fft_chitemp = np.fft.fftn(chitemp)

        mult_result = fft_chitemp * self.dipole_kernel # Elementwise multiplication in freq. domain
        field = np.real(np.fft.ifftn(mult_result))

        field = field[:self.dimensions[0], :self.dimensions[1], :self.dimensions[2]]

        del chitemp


        # To simulate the data acquisition we need to use the signal equation
        # Given we have different labels, the M0 R1 and R2 values can be used from literature

        # We want to create a new volume
        vol = np.zeros(self.dimensions)






    #This is a method that can be later tested but as of now, It should work just fine
    def compute_Bz(self, res, buffer):
        #For later implementation

        # res is resolution from the img header => pixdim
        # Buffer is an int to rescale kspace

        matrix = self.dimensions
   
        # Creating k-space grid
        dim = buffer*matrix

        kmax = 1/(2*res)
        interval = 2*kmax/dim

        kx_ = np.arange(-kmax[0], kmax[0], interval[0])
        ky_ = np.arange(-kmax[1], kmax[1], interval[1])
        kz_ = np.arange(-kmax[2], kmax[2], interval[2])

        kx,ky,kz = np.meshgrid(kx_, ky_, kz_)

        # FFT kernel
        k2 = kx**2 + ky**2 + kz**2
        kernel = np.fft.fftshift(1/3 - kz**2/k2)
        kernel[0,0,0] = 1/3

        FT_chi = np.fft.fftn(self.sus_dist, dim)
        Bz_fft = kernel*FT_chi

        # retrive the inital FOV

        volume_buffed = np.real(np.fft.ifftn(Bz_fft))
        volume_without_buff = volume_buffed[0:matrix[0], 0:matrix[1], 0:matrix[2]]


    def save_sus_csv(self):
        data = []
        for i in self.segmentation_labels.keys():
            label = self.segmentation_labels[i]
            if label.name is not None and label.susceptibility is not None and label.name not in data[1]:
                # The last is to get unique names 
                data.append({"Label ID": label.label_id,
                    "Name": label.name,
                    "Susceptibility": label.susceptibility})
        # Call funtion that creates CSV
        to_csv_sus(data,"data/susceptibility_values.csv")

    def save_relax_csv(self):
        # Further implementation to go through self.relax values of each label?
        # Think about a more efficient way because the user should be able to change the values
        # It might be usefull to get this inputs from different researchers and testing
        pass

    def __repr__(self):
        return f"SegmentationLabelManager == Volume"
