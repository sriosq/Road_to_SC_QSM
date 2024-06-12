import nibabel as nib
import os
import numpy as np
from simulation_functions import *

from label import SegmentationLabel
from Volume import Volume

def plot_kernel_slices(kernel, title="Dipole Kernel Slices"):
    center_slices = [dim // 2 for dim in kernel.shape]

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    axes[0].imshow(np.abs(kernel[center_slices[0], :, :]), cmap='gray')
    axes[0].set_title('Slice at center x')

    axes[1].imshow(np.abs(kernel[:, center_slices[1], :]), cmap='gray')
    axes[1].set_title('Slice at center y')

    axes[2].imshow(np.abs(kernel[:, :, center_slices[2]]), cmap='gray')
    axes[2].set_title('Slice at center z')

    for ax in axes:
        ax.axis('off')

    plt.suptitle(title)
    plt.show()

def create_dipole_kernel2(B0_dir, voxel_size, N, createInKspace):
    B0_dir = np.array(B0_dir)
    B0_dir = B0_dir / np.linalg.norm(B0_dir)
    Nx, Ny, Nz = N
    if createInKspace=='Yes':
        print('Dipole kernel created directly in k-space - it does support rotations')
        kx = np.fft.fftfreq(Nx, voxel_size[0])
        ky = np.fft.fftfreq(Ny, voxel_size[1])
        kz = np.fft.fftfreq(Nz, voxel_size[2])
        kx, ky, kz = np.meshgrid(kx, ky, kz, indexing='ij')
        k_dot_B0 = kx * B0_dir[0] + ky * B0_dir[1] + kz * B0_dir[2]
        k_squared = kx**2 + ky**2 + kz**2
        k_squared[0, 0, 0] = 1  # Avoid division by zero at the origin
        D = 1 / 3 - (k_dot_B0**2 / k_squared)
        D[0, 0, 0] = 0  # Set the DC component to zero
        D = np.fft.fftshift(D)
    else:
        print('Dipole kernel created in the object domain')
        Y, X, Z = np.meshgrid(np.arange(-Ny//2, Ny//2),
                              np.arange(-Nx//2, Nx//2),
                              np.arange(-Nz//2, Nz//2), indexing='ij')
        X = X * voxel_size[0]
        Y = Y * voxel_size[1]
        Z = Z * voxel_size[2]
        d = np.prod(voxel_size) * (3 * (X * B0_dir[0] + Y * B0_dir[1] + Z * B0_dir[2])**2 - X**2 - Y**2 - Z**2) / (4 * np.pi * (X**2 + Y**2 + Z**2)**2.5)
        d[np.isnan(d)] = 0
        D = np.fft.fftn(np.fft.fftshift(d))
    return D

np.set_printoptions(precision=2, suppress=True)


ct_wb_seg = nib.load("post_processing/final_total_seg.nii.gz")

# These variables are instances of the nibabel image
# To work with the data:
# ct_wb_data = ct_wb.get_fdata()
labels = ct_wb_seg.get_fdata()

## Step 1: Creating the dipole kernel
#print(ct_wb_seg.header['pixdim'][1:4])

B0_dir = [0, 0, 1]
voxel_size = ct_wb_seg.header['pixdim'][1:4]
dim = np.array(labels.shape)
dims = [2 * d for d in dim]

D = create_dipole_kernel2(B0_dir, voxel_size, dims,"no")

plot_kernel_slices(D, "Kernel in Object domain")
print(D.shape)

#print("Dipole kernel shape:", D.shape)

# Step 1: complete

#print(ct_wb_seg.header)



