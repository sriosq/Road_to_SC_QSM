# This python file contains functions needed to simulte GRE acquistion

# Creating a dipole kernel
import numpy as np
import matplotlib.pyplot as plt
import nibabel as nib

def create_dipole_kernel(B0_dir, voxel_size,dimensions,buffer=1):

    #B0_dir = np.array(B0_dir)
    #B0_dir = B0_dir / np.linalg.norm(B0_dir)
    dims = dimensions*buffer
    # Extract dimensions for each axis
    kmax = 1/(2*voxel_size)
    interval = 2*kmax/dims
    [kx, ky, kz] = np.meshgrid(np.arange(-kmax[0], kmax[0], interval[0]),
                                np.arange(-kmax[1], kmax[1], interval[1]),
                                np.arange(-kmax[2], kmax[2], interval[2]))
    k_sqrd = kx ** 2 + ky ** 2 + kz ** 2

    with np.errstate(divide='ignore', invalid='ignore'):
        kernel = np.fft.fftshift(1/3 - kz**2/k_sqrd)
        kernel[0,0,0] = 1/3

    return kernel

def show_slices(slices):
   """ Function to display row of image slices """
   fig, axes = plt.subplots(1, len(slices))
   for i, slice in enumerate(slices):
       axes[i].imshow(slice.T, cmap="gray", origin="lower")


def generate_signal(pd, T2star, FA, TE, deltaB0, gamma, handedness):
    # This code is for creating the volume dimensions with new dimension for each TE
    if handedness == 'left':
        sign = -1
    elif handedness == 'right':
        sign = 1
    else:
        raise ValueError("Invalid handedness value")
    # Some regions of the volume have 0 as T2 star values so:
    #if T2star == 0:
    #    signal = pd * np.sin(np.deg2rad(FA))
    #else:
    signal = pd * np.sin(np.deg2rad(FA)) * np.exp(-TE / T2star - sign * 1j * gamma * deltaB0 * TE)
    mag = np.abs(signal)
    phase = np.angle(signal)

    return mag,phase