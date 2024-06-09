# This python file contains functions needed to simulte GRE acquistion

# Creating a dipole kernel
import numpy as np

def create_dipole_kernel(B0_dir, voxel_size,dimensions):

    B0_dir = np.array(B0_dir)
    B0_dir = B0_dir / np.linalg(B0_dir)

    # Extract dimensions for each axis
    Nx, Ny, Nz = dimensions

    # Generate the frequenc grid
    kx = np.fft.fftfreq(Nx, voxel_size[0])
    ky = np.fft.fftfreq(Ny, voxel_size[1])
    kz = np.fft.fftfreq(Nz, voxel_size[2])
    kx, ky, kz = np.meshgrid(kx, ky, kz, indexing='ij')

    k_dot_B0 = kx * B0_dir[0] + ky * B0_dir[1] + kz * B0_dir[2]

    k_squared = kx ** 2 + ky ** 2 + kz ** 2

    k_squared[0, 0, 0] = 1 # This to avoid division by zero at the origin

    dipole_kernel = 1 / 3 - (k_dot_B0 ** 2 / k_squared)
    dipole_kernel[0, 0, 0] = 0  # Set the DC component to zero

    return dipole_kernel

