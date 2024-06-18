import numpy as np
from matplotlib import pyplot as plt
from scipy.ndimage import rotate

class Spherical:
    def __init__(self, matrix, image_res, R, sus_diff):
        self.matrix = matrix
        self.image_res = image_res
        self.R  = R
        self.sus_diff = sus_diff

    def mask(self):
        x_ = np.linspace(-(self.matrix[0]-1)/2, (self.matrix[0]-1)/2, self.matrix[0])
        y_ = np.linspace(-(self.matrix[1]-1)/2, (self.matrix[1]-1)/2, self.matrix[1])
        z_ = np.linspace(-(self.matrix[2]-1)/2, (self.matrix[2]-1)/2, self.matrix[2])

        x,y,z = np.meshgrid(x_, y_, z_)

        r = np.sqrt(x**2 + y**2 + z**2)

        return r**2 < self.R**2
    
    def volume(self):
        return np.where(self.mask() == True, self.sus_diff, 0)
        

class Cylindrical:
    def __init__(self, matrix, image_res, R, theta, sus_diff):
        self.matrix = matrix
        self.image_res = image_res
        self.R  = R
        self.sus_diff = sus_diff
        self.theta = theta # rotation angle about the y-axis

    def mask(self):
        x_ = np.linspace(-(self.matrix[0]-1)/2, (self.matrix[0]-1)/2, self.matrix[0])
        y_ = np.linspace(-(self.matrix[1]-1)/2, (self.matrix[1]-1)/2, self.matrix[1])
        z_ = np.linspace(-(self.matrix[2]-1)/2, (self.matrix[2]-1)/2, self.matrix[2])

        x,y,z = np.meshgrid(x_, y_, z_)

        r = x**2 + y**2 

        mask = r <= self.R**2

        # Rotate the cylinder
        return rotate(mask, self.theta*180/np.pi, axes=(0, 2), reshape=False, order=1)
    
    def volume(self):
        return np.where(self.mask() == True, self.sus_diff, 0)

def computed_Bz(matrix, image_res, sus_dist, buffer):
    #Buffer we take only a limited space in k-space 

    # creating the k-space grid
    dim = buffer*matrix

    kmax = 1/(2*image_res)
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

    return volume_without_buff


