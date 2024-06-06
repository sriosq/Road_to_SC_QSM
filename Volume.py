#Dependencies
import numpy as np
from label import SegmentationLabel

class Volume:
    #establcer relacion volumen label
    
    def __init__(self, volume):
        # we will recieve a 3D so we will follow the convention of labeling for full body scans
        # The 3D volume data will be a np array? or maybe it should be a nifti file and then I can manually
        # get_fdata?

        self.volume = volume
        self.dimensions = np.array(volume.shape)
        self.uniq_labels = np.unique(volume)
        self.segmentation_labels = {} 
        self.sus_dist = np.zeros(self.dimensions)
        # The dictionary has keys for every id number and each value 
        # is the corresponding SegmentationLabel daughter class

        self.create_segmentation_labels()

    #def set_label_name(self, label_id, name):
        #for label in self.segmentation_labels:
            #if label.label_id == label_id:
                #label.set_name(name)
                
    def create_segmentation_labels(self):
        for label_id in self.uniq_labels:
            self.segmentation_labels[label_id] = SegmentationLabel(label_id)

    def set_label_name(self, label_id, name):
        if label_id in self.segmentation_labels:
            self.segmentation_labels[label_id].set_name(name)
        else:
            print(f"Label ID {label_id} not found.")

    def set_label_susceptibility(self, label_id, susceptibility):
        if label_id in self.segmentation_labels:
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


    def compute_Bz(self, res, buffer):
        # res is resolution from the img heaer => pixdim
        # Buffer is an int to rescale kspace

        matrix =np.array(self.volume.shape) # It is initally a tuple but it needs to be an array
   
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

        FT_chi = np.fft.fftn(sus_dist, dim)
        Bz_fft = kernel*FT_chi

        # retrive the inital FOV

        volume_buffed = np.real(np.fft.ifftn(Bz_fft))
        volume_without_buff = volume_buffed[0:matrix[0], 0:matrix[1], 0:matrix[2]]




