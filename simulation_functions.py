# This python file contains functions needed to simulte GRE acquistion

# Creating a dipole kernel
import numpy as np
import matplotlib.pyplot as plt
import nibabel as nib

def create_dipole_kernel(B0_dir, voxel_size,dimensions):

    B0_dir = np.array(B0_dir)
    B0_dir = B0_dir / np.linalg.norm(B0_dir)

    # Extract dimensions for each axis
    Nx, Ny, Nz = dimensions

    # Generate the frequenc grid
    kx = np.fft.fftfreq(Nx, voxel_size[0]).astype(np.float32)
    ky = np.fft.fftfreq(Ny, voxel_size[1]).astype(np.float32)
    kz = np.fft.fftfreq(Nz, voxel_size[2]).astype(np.float32)

    kx, ky, kz = np.meshgrid(kx, ky, kz, indexing='ij')

    k_dot_B0 = kx * B0_dir[0] + ky * B0_dir[1] + kz * B0_dir[2]

    k_squared = kx ** 2 + ky ** 2 + kz ** 2

    k_squared[0, 0, 0] = 1 # This to avoid division by zero at the origin

    dipole_kernel = 1 / 3 - (k_dot_B0 ** 2 / k_squared)
    dipole_kernel[0, 0, 0] = 0  # Set the DC component to zero

    return dipole_kernel

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
        signal = pd * np.sin(np.deg2rad(FA))
    #else:
    signal = pd * np.sin(np.deg2rad(FA)) * np.exp(-TE / T2star - sign * 1j * gamma * deltaB0 * TE)
    mag = np.abs(signal)
    phase = np.angle(signal)

    return mag,phase