#Dependencies
import numpy as np
from label import SegmentationLabel
import nibabel as nib

class Volume:
    #establcer relacion volumen label
    
    def __init__(self, volume):
        # In this version we correct that the output should be the nifti image
        # This way we can attribute the information from nifti files to the class

        self.nifti = volume # This is now directing to a Nifti file
        self.volume = self.nifti.get_fdata()
        self.dimensions = np.array(self.volume.shape)
        self.uniq_labels = np.unique(self.volume)
        self.segmentation_labels = {} 
        self.sus_dist = np.zeros(self.dimensions)
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

        # In TotalSegmentator this is labeld Spinal Cord but
        # it really is the Spinal Canal = Spinal Cord + CSF
        self.set_label_name(76,"SpinalCanal")
        # For SC is (GM + WM)/2 = -9.055
        # From Eva's code GM = -9.03 and WM = -9.08
        self.set_label_susceptibility(76,-9.055)

        # For the lungs we have [9,10,11,12,13]
        for i in [9,10,11,12,13]:
            self.set_label_name(i,'Lungs')
            self.set_label_susceptibility(i,0.35)

        # For the bones we have a lot more labels
        # Vertebrae and ribs. But all of them can have the same value: -11.1
        # Vertebrae list goes from Sacrum to C1
        vertebra_list = np.arange(22,48)
        rib_list = np.append(np.arange(89,114),[68,69,70,71,74,75,88])
        for i in vertebra_list:
            self.set_label_name(i,"Vertebrae")
            self.set_label_susceptibility(i,-11.1)

        for i in rib_list:
            self.set_label_name(i, "Bone")
            self.set_label_susceptibility(i, -11.1)

        # Last but not least is to give susceptibility values to the organs
        # and soft tissue => susceptibility value of water = -9.05

        sus_water_list = [14, 15, 8, 7, 16]

        self.set_label_name(14, "Esophagus")
        self.set_label_name(15, "Trachea")
        self.set_label_name(7, "AdrenalGland")
        self.set_label_name(8, "AdrenalGland")
        self.set_label_name(16, "Thyroid Gland")

        # Susceptibility value of fat = -8.39

        for i in sus_water_list:
            self.set_label_susceptibility(i, -9.05)

        # Soft tissue == water
        # Inside of body == fat

        self.set_label_name(87, "Brain")

        self.set_label_name(0,'Air') # Outside of brain
        self.set_label_susceptibility(0,0.35)

        # If label has not been set it can be considered as fat
        # susceptibility of fat = -8.39

        self.set_label_name(264,"Fat")
        self.set_label_susceptibility(264,-8.39)





    def set_label_name(self, label_id, name):

        if label_id in self.uniq_labels:


            self.segmentation_labels[label_id].set_name(name)
        else:
            print(f"Label ID {label_id} not found.")

    def set_label_susceptibility(self, label_id, susceptibility):
        if label_id in self.uniq_labels:
            # SImilar to set_label_name
            self.segmentation_labels[label_id].set_susceptibility(susceptibility)
        else:
            print(f"Label ID {label_id} not found.")
 
    def display(self):
        pass


    def manual_labeling(self):
        #To iterate over self unique labels we need its total length
        for i in range(len(self.uniq_labels)):
            name = input(f"Enter name for label #{i}: ") 
            label = SegmentationLabel(i, name)



    def show_labels(self):
        for i in range(len(self.uniq_labels)):
            print(i)

    def __repr__(self):
        return f"SegmentationLabelManager == Volume"

    def create_sus_dist(self):
        
        dimensions =np.array(self.volume.shape)

        for i in range(dimensions[0]):
            for j in range(dimensions[1]):
                for k in range(dimensions[2]):

                    pixel = self.volume[i,j,k]
                    label = self.segmentation_labels[pixel]
                    suscep = label.susceptibility
                    self.sus_dist[i,j,k] = suscep

        return self.sus_dist

    def save_sus_dist_nii(self):
        # Method to save the susceptibility distribution created to nifti
        new_img = nib.Nifti1Image(self.sus_dist, affine=self.nifti.affine)

        # Save the new NIfTI image to a file

        nib.save(new_img,"sus_dist.nii.gz")

    def compute_Bz(self, res, buffer):
        #For later implementation

        # res is resolution from the img header => pixdim
        # Buffer is an int to rescale kspace

        matrix =np.array(self.volume.shape) # It is initially a tuple, but it needs to be an array
   
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




