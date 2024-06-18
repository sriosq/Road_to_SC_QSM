#Dependencies
import numpy as np
from label import SegmentationLabel
import nibabel as nib
from simulate_mri import *
from utils.simulation_functions import *
import scipy.ndimage
from utils.get_dic_values import *

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
        self.pd_dist = np.zeros(self.dimensions)
        self.deltaB0 = np.zeros(self.dimensions)

        self.magnitude = np.zeros(self.dimensions)
        self.phase = np.zeros(self.dimensions)
        # The dictionary has keys for every id number and each value 
        # is the corresponding SegmentationLabel daughter class

        self.create_labels()
        # This is the Convention Dictionary for labels_id - names - sus_values

    def create_labels(self):
       for label_id in self.uniq_labels:
           self.segmentation_labels[label_id] = SegmentationLabel(label_id)

    def create_segmentation_labels(self):

        # The most important labels are:
        # lungs Bone Soft Tissue SpinalCord CSF
        # These are the regions of the body that impact the most
        #
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
            self.set_label_susceptibility(i,0.35)

        # For the bones we have a lot more labels
        # Vertebrae and ribs. But all of them can have the same value: -11.1
        # Vertebrae list goes from Sacrum to C1
        vertebra_list = np.arange(22,48)
        rib_list = np.append(np.arange(89,114),[68,69,70,71,74,75,88])
        for i in vertebra_list:
            self.set_label_name(i,"bone")
            self.set_label_susceptibility(i,-11.1)

        for i in rib_list:
            self.set_label_name(i, "bone")
            self.set_label_susceptibility(i, -11.1)

        # Last but not least is to give susceptibility values to the organs
        # and soft tissue => susceptibility value of water = -9.05

        sus_water_list = [1, 2, 3, 4, 5, 6, 7, 8, 14, 15, 16]

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

        # Soft tissue == water
        # Inside of body == fat


        self.set_label_name(0,"air") # Outside of brain
        self.set_label_susceptibility(0,0.35)

        # If label has not been set it can be considered as fat
        # susceptibility of fat label = (sus_fat + sus_muscle)/2 = -8.39 + -9.03 div2 = -8.71

        self.set_label_name(264,"fat")
        self.set_label_susceptibility(264,-8.71)

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
        return self.pd_dist
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
        pass


    def manual_labeling(self):
        #To iterate over self unique labels we need its total length
        for i in range(len(self.uniq_labels)):
            name = input(f"Enter name for label #{i}: ") 
            label = SegmentationLabel(i, name)
            # Complete later for all attributes

    def show_labels(self):
        for i in range(len(self.uniq_labels)):
            print(i)


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
        # Save the new NIfTI image to a file
        nib.save(temp_img,"sus_dist.nii.gz")
        del temp_img

    def save_pd_dist(self):
        # Method to save the proton density distribution created to nifti
        temp_img = nib.Nifti1Image(self.pd_dist, affine=self.nifti.affine)
        # Save the new NIfTI image to a file
        nib.save(temp_img,"pd_dist.nii.gz")
        del temp_img

    def calculate_deltaB0(self,B0_dir =[0,0,1]):

        voxel_size = self.nifti.header["pixdim"][1:4]
        padded_dims = tuple(2*self.dimensions)
        D = create_dipole_kernel(B0_dir,voxel_size,self.dimensions)

        sus_dist_padded = np.zeros(padded_dims,dtype=np.float32)
        sus_dist_padded[:self.dimensions[0], :self.dimensions[1], :self.dimensions[2]] = self.sus_dist

        self.deltaB0 = np.real(np.fft.ifftn(np.fft.fftn(sus_dist_padded)*D))

    def save_deltaB0(self):
        temp_img = nib.Nifti1Image(self.deltaB0, affine=self.nifti.affine)
        nib.save(temp_img,"dipole_kernel.nii.gz")
        del temp_img
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
        self.measurement = np.zeros(newVol_dims)
        gamma = 267.52218744e6*3/B0# Using gamma for 3 Tesla, B0 is optional to change
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

        return self.measurement
    # Line to implement from MATLAB
    # signal = pd*sind(FA)*exp(-TE./T2star-sign*1i*gamma*deltaB0*TE)

    def get_Magnitude(self):
        self.magnitude = np.abs(self.measurement)
        temp_img = nib.Nifti1Image(self.magnitude, affine=self.nifti.affine)
        # Save the new NIfTI image to a file
        nib.save(temp_img, "magnitude_simulated.nii.gz")
        del temp_img

    def get_Phase(self):
        self.phase = np.angle(self.measurement)
        temp_img = nib.Nifti1Image(self.phase, affine=self.nifti.affine)
        # Save the new NIfTI image to a file
        nib.save(temp_img, "phase_simulated.nii.gz")
        del temp_img

    def get_Real(self):
        self.real = np.real(self.measurement)
        temp_img = nib.Nifti1Image(self.real, affine=self.nifti.affine)
        # Save the new NIfTI image to a file
        nib.save(temp_img, "real_simulated.nii.gz")
        del temp_img

    def get_Imaginary(self):
        self.imag = np.imag(self.measurement)
        temp_img = nib.Nifti1Image(self.imag, affine=self.nifti.affine)
        # Save the new NIfTI image to a file
        nib.save(temp_img, "imaginary_simulated.nii.gz")
        del temp_img



    # Implementation of the code to simulate MRI data acquisition
    # This code is more complex as it assumes that there is still Longitudinal Magnetizatation
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
        to_csv_sus(data,"susceptibility_values.csv")

    def save_relax_csv(self):
        # Further implementation to go through self.relax values of each label?
        # Think about a more efficient way because the user should be able to change the values
        # It might be usefull to get this inputs from different researchers and testing
        pass

    def __repr__(self):
        return f"SegmentationLabelManager == Volume"
